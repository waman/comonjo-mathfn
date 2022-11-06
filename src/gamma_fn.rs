use crate::integer_util::*;

const PI: f64 = std::f64::consts::PI;
const LOG_2PI: f64 = 1.8378770664093453_f64;

const N: f64 = 8.;

const B2 : f64 = 1. / 6.;
const B4 : f64 = -1. / 30.;
const B6 : f64 = 1. / 42.;
const B8 : f64 = -1. / 30.;
const B10: f64 = 5. / 66.;
const B12: f64 = -691. / 2730.;
const B14: f64 = 7. / 6.;
const B16: f64 = -3617. / 510.;

/// Return the logarithm of the gamma function *log Γ(x)*.
/// (The argument *x* must be positive.)
///
/// Ref: 『改訂新版 Cによる標準アルゴリズム事典』ガンマ関数 (gamma function) gamma.c
pub fn log_gamma(mut x: f64) -> f64 {
    if x.is_infinite() { return if x > 0. { f64::INFINITY }else{ f64::NAN }; }

    let mut v: f64 = 1.;
    while x < N {
        v *= x;
        x += 1.;
    }

    let w = 1. / (x * x);
    ((((((((B16 / (16. * 15.))  * w + (B14 / (14. * 13.))) * w
         + (B12 / (12. * 11.))) * w + (B10 / (10. *  9.))) * w
         + (B8  / ( 8. *  7.))) * w + (B6  / ( 6. *  5.))) * w
         + (B4  / ( 4. *  3.))) * w + (B2  / ( 2. *  1.))) / x
         + 0.5 * LOG_2PI - v.ln() - x + (x - 0.5) * x.ln()
}

#[test]
fn test_the_values_of_log_gamma_at_the_specific_points(){
    
    assert!(log_gamma(1.).abs() < 1e-12, "logΓ(1) = 0");
    
    // Non-positive integers
    assert!(log_gamma(0.).is_infinite(),  "logΓ(0) = ∞");
    assert!(log_gamma(-1.).is_infinite(), "logΓ(-1) = ∞");
    assert!(log_gamma(-2.).is_infinite(), "logΓ(-2) = ∞");
    
    // Negative non-integer numbers
    assert!(log_gamma(-0.5).is_nan(), "logΓ(-0.5) = NaN");
    assert_eq!(log_gamma(-1.5), 0.8600470153764803_f64, "logΓ(-1.5) = logΓ(0.5) - ln(1.5) - ln(0.5)");
    assert!(log_gamma(-2.5).is_nan(), "logΓ(-2.5) = NaN");
}

#[test]
fn test_the_values_of_log_gamma_at_non_finite_args(){
    assert_eq!(log_gamma(f64::INFINITY), f64::INFINITY, "logΓ(∞) = ∞");
    assert!(log_gamma(f64::NEG_INFINITY).is_nan(),      "logΓ(-∞) = NaN");
    assert!(log_gamma(f64::NAN).is_nan(),               "logΓ(NaN) = NaN");
}

/// Return a value of the gamma function *Γ(x)*.
/// 
/// Ref: 『改訂新版 Cによる標準アルゴリズム事典』ガンマ関数 (gamma function) gamma.c
/// 
/// #### The value at the special points ####
/// 
/// ```
/// use comonjo_mathfn::gamma;
/// 
/// fn assert_approximately(x: f64, y: f64){
///     assert!((x - y).abs() <= 1e-14);
/// }
/// 
/// assert_approximately(gamma(1.), 1.);
/// assert_approximately(gamma(2.), 1.);
/// assert_approximately(gamma(3.), 2. * 1.);
/// assert_approximately(gamma(4.), 3. * 2. * 1.);
/// assert_approximately(gamma(5.), 4. * 3. * 2. * 1.);
/// 
/// use f64::consts;
/// assert_approximately(gamma(0.5), consts::PI.sqrt());
/// ```
/// 
/// #### Diverging points ####
/// ```
/// use comonjo_mathfn::gamma;
/// 
/// assert!(gamma(0.).is_nan());
/// assert!(gamma(-1.).is_nan());
/// assert!(gamma(-2.).is_nan());
/// assert!(gamma(-3.).is_nan());
/// assert!(gamma(-4.).is_nan());
/// assert!(gamma(-5.).is_nan());
/// 
/// //assert!(gamma(0.).is_infinite());
/// //assert!(gamma(-1.).abs() > 1e16);
/// //assert!(gamma(-2.).abs() > 1e15);
/// //assert!(gamma(-3.).abs() > 1e15);
/// //assert!(gamma(-4.).abs() > 1e14);
/// //assert!(gamma(-5.).abs() > 1e13);
/// ```
pub fn gamma(x: f64) -> f64 {
    if x > 0. {
        log_gamma(x).exp()
    }else if x.is_integer() {
        f64::NAN
    }else{
        PI / ((PI * x).sin() * log_gamma(1. - x).exp())
    }
}

#[cfg(test)]
use super::test_util::*;

#[test]
fn test_const_value(){
    assert_eq!(LOG_2PI, std::f64::consts::TAU.ln());
}

#[test]
fn test_the_values_of_gamma_at_positive_integers_are_factorials(){
    let mut exp = 1.;
    for n in 1..20 {
        let x = n as f64;
        assert_approximately(gamma(x), exp, EPS, "Γ(n) = (n-1)!");
        exp *= x;
    }
}

#[test]
fn test_the_values_of_gamma_at_half_integers(){
    let sqrt_pi: f64 = PI.sqrt();
    assert_approximately(gamma(0.5), sqrt_pi,       EPS, "Γ(1/2) = √π");
    assert_approximately(gamma(1.5), sqrt_pi/2.,    EPS, "Γ(3/2) = √π/2");
    assert_approximately(gamma(2.5), sqrt_pi*3./4., EPS, "Γ(5/2) = 3√π4");
    
    assert_approximately(gamma(-0.5), -2.*sqrt_pi,   EPS, "Γ(-1/2) = -2√π");
    assert_approximately(gamma(-1.5), 4.*sqrt_pi/3., EPS, "Γ(-3/2) = 4√π/3");
}

#[test]
fn test_that_the_gamma_function_diverges_at_zero_or_negative_integers(){
    for n in 0..=5 {
    let x = -n as f64;
        assert_diverging(gamma(x), BIG_VALUE, &format!("Γ(-n) ~ ±∞ at -n = {}", x));
    }
}

#[test]
fn test_the_values_of_gamma_at_non_finite_args(){
    assert_eq!(gamma(f64::INFINITY), f64::INFINITY, "Γ(∞) = ∞");
    assert!(gamma(f64::NEG_INFINITY).is_nan(),      "Γ(-∞) = NaN");
    assert!(gamma(f64::NAN).is_nan(),               "Γ(NaN) = NaN");
}

#[test]
fn test_the_gamma_function_properties(){
    const DELTA: f64 = 1e-3;

    should_the_same_mathfn(
        "Γ(x + 1) = xΓ(x)", 
            |x| gamma(x + 1.), 
            |x| x * gamma(x))
        .filter(|x| !is_close_to_a_non_positive_integer(x, DELTA)).assert();

    should_the_same_mathfn(
        "Γ(1 - x)Γ(x) = π/sin(πx)", 
            |x| gamma(1. - x) * gamma(x), 
            |x| PI / (PI * x).sin())
        .filter(|x| !is_close_to_an_integer(x, DELTA))
        .epsilon(1e-11).assert();
}