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
/// assert_approximately(beta(0.5, 0.5), std::f64::consts::PI);
/// 
/// assert_approximately(beta(1., 1.), 1.);
/// 
/// assert_approximately(beta(2., 1.), 0.5);
/// assert_approximately(beta(1., 2.), 0.5);
/// 
/// assert_approximately(beta(3., 1.), 1./3.);
/// assert_approximately(beta(2., 2.), 1./6.);
/// assert_approximately(beta(1., 3.), 1./3.);
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
        (log_gamma(x) + log_gamma(y) - log_gamma(x + y)).exp()
    }else{
        gamma(x) * gamma(y) / gamma(x + y)
    }
}

#[cfg(test)]
use crate::test_util::*;
#[cfg(test)]
use std::f64::consts::PI;

#[test]
fn test_diverging_if_the_one_of_the_arguments_is_zero_or_negative_integers(){
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
fn test_the_beta_function_properties(){
    const DELTA: f64 = 1e-3;
    
    should_the_same_mathfn2(
        "Β(x, y) = Β(y, x)", 
            |x, y| beta(x, y), 
            |x, y| beta(y, x))
        .filter(|x, y| !x.is_non_positive_integer(DELTA) && !y.is_non_positive_integer(DELTA))
        .assert();

    should_the_same_mathfn2(
        "xΒ(x, y+1) = yΒ(x+1, y)", 
            |x, y| x * beta(x, y + 1.), 
            |x, y| y * beta(x + 1., y))
        .filter(|x, y| !x.is_non_positive_integer(DELTA) && !y.is_non_positive_integer(DELTA))
        .epsilon(1e-11).assert();

    should_the_same_mathfn2(
        "B(x, y) = B(x, y+1) + B(x+1, y)", 
            |x, y| beta(x, y), 
            |x, y| beta(x, y + 1.) + beta(x + 1., y))
        .filter(|x, y| !x.is_non_positive_integer(DELTA) && !y.is_non_positive_integer(DELTA))
        .epsilon(1e-11).assert();

    should_the_same_mathfn(
        "B(1, x) = 1/x", 
            |x| beta(1., x), 
            |x| 1. / x)
        .filter(|x| x > 0.).assert();

    should_the_same_mathfn(
        "B(x, 1-x) = π/sin(πx)", 
            |x| beta(x, 1. - x), 
            |x| PI / (PI * x).sin())
        .filter(|x| !x.is_integer(DELTA)).assert();
}
