use crate::gamma_fn::{log_gamma, gamma};

/// Return a value of the beta function *Β(x, y)*.
/// 
/// Ref: 『改訂新版 Cによる標準アルゴリズム事典』ベータ関数 (beta function) beta.c
/// 
/// #### The value at the special points ####
/// 
/// ```
/// use comonjo_mathfn::beta;
/// 
/// fn assert_approximately(x: f64, y: f64){
///     assert!((x - y).abs() <= 1e-14);
/// }
/// 
/// assert_approximately(beta(0.5, 0.5), f64::consts::PI);
/// 
/// assert_approximately(beta(1., 1.), 1.);
/// 
/// assert_approximately(beta(2., 1.), 0.5);
/// assert_approximately(beta(1., 2.), 0.5);
/// 
/// assert_approximately(beta(3., 1.), 1./3.);
/// assert_approximately(beta(2., 2.), 1./6.);
/// assert_approximately(beta(1., 3.), 1./3.);
/// 
/// assert_approximately(beta(-3., 1.), -2.);
/// ```
/// 
/// #### Diverging points ####
/// ```
/// use comonjo_mathfn::beta;
/// 
/// assert!(beta(0., 0.).is_nan());
/// assert!(beta(1., 0.).is_infinite());
/// assert!(beta(0., 1.).is_infinite());
/// ```
pub fn beta(x: f64, y: f64) -> f64 {
    if x > 0. && y > 0. {
        if x.is_infinite() || y.is_infinite() {
            0.
        }else{
            (log_gamma(x) + log_gamma(y) - log_gamma(x + y)).exp()
        }
    }else{
        if x.is_infinite() || y.is_infinite() {
            if x.is_infinite() && y.is_infinite() {
                // never occurs x > 0 && y > 0
                if x < 0. && y < 0. { f64::INFINITY }else{ f64::NAN }
            }else{
                let inf = if x.is_infinite() { x }else{ y };
                if inf > 0. { 0. }else{ f64::INFINITY }
            }
        }else{
            gamma(x) * gamma(y) / gamma(x + y)
        }
    }
}

#[cfg(test)]
use crate::test_util::*;
#[cfg(test)]
use std::f64::consts::PI;

#[test]
fn test_that_the_beta_function_diverges_if_the_one_of_the_arguments_is_zero_or_negative_integers(){
    repeat(10, ||{
        let y: f64 = rand(-5., 5.);
        let big_value: f64 = if (y - y.round()).abs() < 0.05 { 1e8 } else { BIG_VALUE };
        for i in 0..=5 {
            let x = -i as f64;

            assert_diverging(beta(x, y), big_value,
                &format!("B(-n, y) ~ ±∞ at (-n, y) = ({}, {})", x, y));

            assert_diverging(beta(y, x), big_value,
                &format!("B(x, -n) ~ ±∞ at (x, -n) = ({}, {})", y, x));
        }
    })
}

#[test]
fn test_the_values_of_beta_at_the_infinities_and_nan(){

    let values = &[f64::NAN, f64:: NEG_INFINITY, f64::INFINITY, 0.5];
    values.iter().zip(values.iter()).for_each(|(ra, rb)|{
        let (a, b) = (*ra, *rb);

        if a.is_nan() || b.is_nan() {
            assert!(beta(a, b).is_nan(), "Β({}, {}) = NaN", a, b);

        }else if a.is_infinite() && b.is_infinite() {
            if a > 0. && b > 0. {
                assert_approximately(beta(a, b), 0., EPS, &format!("Β({}, {}) = 0", a, b));
            }else if a < 0. && b < 0. {
                assert!(beta(a, b) == f64::INFINITY, "Β({}, {}) = ∞", a, b);
            }else{
                assert!(beta(a, b).is_nan(), "Β({}, {}) = NaN", a, b);
            }

        }else if a.is_infinite() || b.is_infinite() {
            let inf = if a.is_infinite() { a } else { b };
            if inf > 0. {
                assert_approximately(beta(a, b), 0., EPS, &format!("Β({}, {}) = 0", a, b));
            }else{
                assert!(beta(a, b) == f64::INFINITY, "Β({}, {}) = ∞", a, b);
            }

        }else{
            assert!(beta(a, b).is_finite(), "|Β({}, {})| < ∞", a, b);
        }
    });
}

#[test]
fn test_the_beta_function_properties(){
    const DELTA: f64 = 1e-3;
    
    should_the_same_mathfn2(
        "Β(x, y) = Β(y, x)", 
            |x, y| beta(x, y), 
            |x, y| beta(y, x))
        .filter(|x, y| !is_close_to_a_non_positive_integer(x, DELTA) && !is_close_to_a_non_positive_integer(y, DELTA))
        .assert();

    should_the_same_mathfn2(
        "xΒ(x, y+1) = yΒ(x+1, y)", 
            |x, y| x * beta(x, y + 1.), 
            |x, y| y * beta(x + 1., y))
        .filter(|x, y| [x, y, x+y].iter().all(|t| !is_close_to_a_non_positive_integer(*t, DELTA)))
        .epsilon(1e-11).assert();

    should_the_same_mathfn2(
        "B(x, y) = B(x, y+1) + B(x+1, y)", 
            |x, y| beta(x, y), 
            |x, y| beta(x, y + 1.) + beta(x + 1., y))
        .filter(|x, y| [x, y, x+y].iter().all(|t| !is_close_to_a_non_positive_integer(*t, DELTA)))
        .epsilon(1e-11).assert();

    should_the_same_mathfn(
        "B(1, x) = 1/x", 
            |x| beta(1., x), 
            |x| if !is_close_to_a_non_positive_integer(x, 1e-5) {
                    1. / x
                }else{ 
                    f64::NAN
                }).assert();

    should_the_same_mathfn(
        "B(x, 1-x) = π/sin(πx)", 
            |x| beta(x, 1. - x), 
            |x| PI / (PI * x).sin())
        .filter(|x| !is_close_to_an_integer(x, DELTA)).assert();
}
