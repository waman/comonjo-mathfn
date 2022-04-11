const PI: f64 = std::f64::consts::PI;
const LOG_2PI: f64 = 1.8378770664093453_f64;

const N: f64 = 8.0;

const B2 : f64 = 1.0 / 6.0;
const B4 : f64 = -1.0 / 30.0;
const B6 : f64 = 1.0 / 42.0;
const B8 : f64 = -1.0 / 30.0;
const B10: f64 = 5.0 / 66.0;
const B12: f64 = -691.0 / 2730.0;
const B14: f64 = 7.0 / 6.0;
const B16: f64 = -3617.0 / 510.0;

/// Return the logarithm of the gamma function *log Γ(x)*.
/// (The argument *x* must be positive.)
///
/// Ref: 『改訂新版 Cによる標準アルゴリズム事典』ガンマ関数 (gamma function) gamma.c
fn log_gamma(mut x: f64) -> f64 {
    debug_assert!(x > 0.0, "argument must be positive: {}", x);

    let mut v: f64 = 1.0;
    while x < N {
        v *= x;
        x += 1.0;
    }

    let w = 1.0 / (x * x);
    ((((((((B16 / (16.0 * 15.0))  * w + (B14 / (14.0 * 13.0))) * w
         + (B12 / (12.0 * 11.0))) * w + (B10 / (10.0 *  9.0))) * w
         + (B8  / ( 8.0 *  7.0))) * w + (B6  / ( 6.0 *  5.0))) * w
         + (B4  / ( 4.0 *  3.0))) * w + (B2  / ( 2.0 *  1.0))) / x
         + 0.5 * LOG_2PI - v.ln() - x + (x - 0.5) * x.ln()
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
/// fn assert_equiv(x: f64, y: f64){
///     assert!((x - y).abs() <= 1e-14);
/// }
/// 
/// assert_equiv(gamma(1.0), 1.0);
/// assert_equiv(gamma(2.0), 1.0);
/// assert_equiv(gamma(3.0), 2.0 * 1.0);
/// assert_equiv(gamma(4.0), 3.0 * 2.0 * 1.0);
/// assert_equiv(gamma(5.0), 4.0 * 3.0 * 2.0 * 1.0);
/// 
/// use std::f64::consts;
/// assert_equiv(gamma(0.5), consts::PI.sqrt());
/// ```
/// 
/// #### Diverging points ####
/// ```
/// use comonjo_mathfn::gamma;
/// 
/// assert!(gamma(0.0).is_infinite());
/// assert!(gamma(-1.0).abs() > 1e16);
/// assert!(gamma(-2.0).abs() > 1e15);
/// assert!(gamma(-3.0).abs() > 1e15);
/// assert!(gamma(-4.0).abs() > 1e14);
/// assert!(gamma(-5.0).abs() > 1e13);
/// ```
pub fn gamma(x: f64) -> f64 {
    if x <= 0.0 {
        PI / ((PI * x).sin() * log_gamma(1.0 - x).exp())
    }else{
        log_gamma(x).exp()
    }
}

#[cfg(test)]
use super::test_util::*;

#[test]
fn test_consts(){
    assert_eq!(LOG_2PI, std::f64::consts::TAU.ln());
}

#[test]
fn test_the_values_of_positive_integers_are_factorials(){
    let mut exp = 1.0;
    for n in 1..20 {
        let x = n as f64;
        assert_approximately(gamma(x), exp, EPS, "Γ(n) = (n-1)!");
        exp *= x;
    }
}

#[test]
fn test_the_values_of_half_integers(){
    let sqrt_pi: f64 = PI.sqrt();
    assert_approximately(gamma(0.5), sqrt_pi,         EPS, "Γ(1/2) = √π");
    assert_approximately(gamma(1.5), sqrt_pi/2.0,     EPS, "Γ(3/2) = √π/2");
    assert_approximately(gamma(2.5), sqrt_pi*3.0/4.0, EPS, "Γ(5/2) = 3√π4");
    
    assert_approximately(gamma(-0.5), -2.0*sqrt_pi,    EPS, "Γ(-1/2) = -2√π");
    assert_approximately(gamma(-1.5), 4.0*sqrt_pi/3.0, EPS, "Γ(-3/2) = 4√π/3");
}

#[test]
fn test_diverging_at_zero_or_negative_integers(){
    for n in 0..=5 {
        let x = -n as f64;
        assert_diverging(gamma(x), BIG_VALUE, "Γ(-n) ~ ±∞");
    }
}

#[test]
fn test_the_gamma_function_properties(){
    assert_the_same_mathfn(
        |x| gamma(x + 1.0), 
        |x| x * gamma(x), 
        |x| !(x <= 0.0 && x.round() == x), "Γ(x + 1) = xΓ(x)");

    assert_the_same_mathfn(
        |x| gamma(1.0 - x) * gamma(x), 
        |x| PI / (PI * x).sin(), 
        |x| x.round() != x, "Γ(1 - x)Γ(x) = π/sin(πx)");
}