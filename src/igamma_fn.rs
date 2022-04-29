// The tests refer to
// <a href="https://en.m.wikipedia.org/wiki/Incomplete_gamma_function">Incomplete gamma function</a>

use crate::gamma_fn::log_gamma;

/** (log π)/2 */
const LOG_PI_BY2: f64 = 0.5723649429247001;

/// Return a value of the the regularized incomplete gamma function.
/// The normalization factor can be manually specified by the last argument.
/// 
/// Ref: 『改訂新版 Cによる標準アルゴリズム事典』不完全ガンマ関数 (incomplete gamma function) igamma.c
fn p_gamma_normalizable(s: f64, x: f64, log_gamma_s: f64) -> f64 {
    if x >= 1. + s { return 1. - q_gamma_normalizable(s, x, log_gamma_s); }
    if x == 0. { return 0.; }

    let mut result = x.powf(s) * (-x - log_gamma_s).exp() / s;
    // let mut result = (s * x.ln() - x - log_gamma_s).exp() / s;
        // The above line doesn't work when s is a non-positive integer and x is negative.
    let mut term = result;
    let mut k = 1.;
    while k < 1000. {
        let prev = result;
        term *= x / (s + k);
        result += term;
        if result == prev { return result; }
        k += 1.;
    }
    // (1..).map(|k| k as f64).map(|k| x / (s + k))
    // .scan(result, |term, factor|{
    //     *term = *term * factor;
    //     Some(*term)
    // }).scan(0., |sum, term|{
    //     *sum = *sum + term;
    //     Some(sum)
    // })
    std::f64::NAN
} 

/// Return a value of the upper incomplete gamma function.
/// The normalization factor can be manually specified by the last argument.
/// 
/// Ref: 『改訂新版 Cによる標準アルゴリズム事典』不完全ガンマ関数 (incomplete gamma function) igamma.c
fn q_gamma_normalizable(s: f64, x: f64, log_gamma_s: f64) -> f64 {
    if x < 1. + s { return 1. - p_gamma_normalizable(s, x, log_gamma_s); }

    let mut w = x.powf(s) * (-x - log_gamma_s).exp();
    // let mut w = (s * x.ln() - x - log_gamma_s).exp();
        // The above line doesn't work when s is a non-positive integer and x is negative.
    let mut la = 1.; let mut lb = 1. + x - s;
    let mut result = w / lb;
    let mut k = 2.;
    while k < 1000. {
        let mut temp = ((k-1.-s)*(lb-la) + (k+x)*lb)/k;
        la = lb; lb = temp;
        w *= (k-1.-s)/k;
        temp = w/(la*lb);
        let prev = result;
        result += temp;
        if result == prev { return result; }
        k += 1.;
    }
    return std::f64::NAN;
}

//***** incomplete gamma function *****
/// Return a value of the lower incomplete gamma function *γ(s, x)*.
pub fn igamma(s: f64, x: f64) -> f64 {
    let log_gamma_s = log_gamma(s);
    log_gamma_s.exp() * p_gamma_normalizable(s, x, log_gamma_s)
}

/// Return a value of the upper incomplete gamma function *Γ(s, x) = Γ(s) - γ(s, x)*.
#[allow(non_snake_case)]
pub fn iGamma(s: f64, x: f64) -> f64 {
    let log_gamma_s = log_gamma(s);
    log_gamma_s.exp() * q_gamma_normalizable(s, x, log_gamma_s)
}

#[test]
fn igamma_converge_when_s_is_a_positive_integer_and_x_is_negative(){
    use crate::test_util::*;
    repeat(10, ||{
        let x = rand(-10., 0.);
        for n in 1..=5 {
            assert!(igamma(n as f64, x).is_finite(), 
                "|γ(n, x)| < ∞ at (n, x) = ({}, {})", n, x);
        }
    });
}

#[test]
fn igamma_diverge_when_s_is_a_non_positive_integer_and_x_is_negative(){
    repeat(10, ||{
        let x = rand(-10., 0.);
        for n in 1..=5 {
            assert_diverging(igamma(-n as f64, x), BIG_VALUE, 
                &format!("γ(-n, x) ~ ±∞ at (-n, x) = ({}, {})", -n, x));
        }
    });
}

#[test]
fn test_the_incomplete_gamma_function_properties(){
    should_the_same_mathfn2(
        "γ(s, x) + Γ(s, x) = Γ(s)",
            |s, x| igamma(s, x) + iGamma(s, x), 
            |s, _| gamma(s)) 
        .filter(|s, x| x < s)
        .var0(|v| v.name("s").range(0., 10.).end())
        .var1(|v| v.name("x").range(0., 10.).end()).assert();

    should_the_same_mathfn2(
        "γ(s+1, x) = s*γ(s, x) - x^s*e^{-x}",
            |s, x| igamma(s+1., x), 
            |s, x| s*igamma(s, x) - x. powf(s) * (-x).exp())
        .filter(|s, x| s != 0. && x != 0.)
        .var0(|v| v.name("s").range(0., 10.).end())
        .var1(|v| v.name("x").range(0., 10.).end()).assert();

    should_the_same_mathfn2(
        "Γ(s+1, x) = s*Γ(s, x) + x^s*e^{-x}",
            |s, x| iGamma(s+1., x), 
            |s, x| s*iGamma(s, x) + x.powf(s) * (-x).exp())
        .filter(|s, x| s != 0. && x != 0.)
        .var0(|v| v.name("s").range(0., 10.).end())
        .var1(|v| v.name("x").range(0., 10.).end()).assert();

    should_the_same_mathfn(
        "γ(1, x) = 1 - e^{-x}",
            |x| igamma(1., x), 
            |x| 1. - (-x).exp()).assert();

    should_the_same_mathfn(
        "Γ(1, x) = e^{-x}", 
            |x| iGamma(1., x), 
            |x| (-x).exp()).assert();

    should_the_same_mathfn(
        "Γ(s, 0) = Γ(s)",
            |s| iGamma(s, 0.), 
            |s| gamma(s))
        .var0(|v| v.name("s").range(0., 20.).end()).assert();

    use std::f64::consts;

    fn factorial(n: f64) -> f64 {
        let mut i = n.round() as usize;
        if i == 1 { return 1.}

        let mut result = 1.;
        while i > 1 {
            result *= i as f64;
            i -= 1;
        }

        result
    }

    should_the_same_mathfn(
        "Γ(s+1, 1) = [e*s!]/e",
            |s| iGamma(s+1., 1.),
            |s| (consts::E * factorial(s)).floor() / consts::E)
        .var0(|v| v.name("s").range(1., 10.).is_integer(true).end()).assert();

    should_the_same_mathfn2(
        "Γ(s, x) = (s-1)!e^{-x}sum_0^{s-1}x^k/k!",
            |s, x| iGamma(s, x),
            |s, x| {
                let sint = s.round() as i32;
                let mut sum = 0.;
                let mut k: i32 = 0;
                let mut kk = 1.;
                while k < sint {
                    if k != 0 { kk *= k as f64 };
                    sum += x.powi(k) / kk;
                    k += 1;
                }
                return factorial(s-1.) * (-x).exp() * sum;
            })
        .var0(|v| v.name("s").range(1., 10.).is_integer(true).end())
        .var1(|v| v.name("x").end())
        .epsilon(1e-10).assert();
}

/// Return a value of the regularized Gamma function *P(s, x) = γ(s, x)/Γ(s)*.
pub fn p_gamma(s: f64, x: f64) -> f64 {
    p_gamma_normalizable(s, x, log_gamma(s))
}

/// Return a value of the regularized Gamma function *Q(s, x)* = Γ(s, x)/Γ(s).
pub fn q_gamma(s: f64, x: f64) -> f64 {
    q_gamma_normalizable(s, x, log_gamma(s))
}

#[test]
fn test_regularized_gamma_function_properties(){

    should_the_same_mathfn2(
        "P(s, x) = γ(s, x)/Γ(s)",
            |s, x| p_gamma(s, x),
            |s, x| igamma(s, x) / gamma(s))
        .filter(|s, _| s > 0.)
        .var0(|v| v.name("s").end())
        .var1(|v| v.name("x").end()).assert();

    should_the_same_mathfn2(
        "Q(s, x) = Γ(s, x)/Γ(s)", 
            |s, x| q_gamma(s, x),
            |s, x| iGamma(s, x) / gamma(s))
        .filter(|s, _| s > 0.)
        .var0(|v| v.name("s").end())
        .var1(|v| v.name("x").end()).assert();
}

//***** Gauss error function *****
/// Return a value of the (Gauss) error function.
/// 
/// Ref: 『改訂新版 Cによる標準アルゴリズム事典』不完全ガンマ関数 (incomplete gamma function) igamma.c
pub fn erf(x: f64) -> f64 {
    if x >= 0. {
        p_gamma_normalizable(0.5, x*x, LOG_PI_BY2)
    }else{
        -p_gamma_normalizable(0.5, x*x, LOG_PI_BY2)
    }
}

/// Return a value of the complementary error function *erfc(x)* (= *1 - erf(x)*).
/// 
/// Ref: 『改訂新版 Cによる標準アルゴリズム事典』不完全ガンマ関数 (incomplete gamma function) igamma.c
pub fn erfc(x: f64) -> f64 {
    if x >= 0. {
        q_gamma_normalizable(0.5, x*x, LOG_PI_BY2)
    }else{
        1. + p_gamma_normalizable(0.5, x*x, LOG_PI_BY2)
    }
}

#[test]
fn test_error_function_properties(){
    let sqrt_pi_inv = 1. / (std::f64::consts::PI.sqrt());

    should_the_same_mathfn(
        "erfc(x) = 1 - erf(x)",
            |x| erfc(x),
            |x| 1. - erf(x))
        .var0(|v| v.range(-3., 3.).end()).assert();

    should_the_same_mathfn(
        "erf(x) = sgn(x)γ(1/2, x^2)/√π",
            |x| erf(x),
            |x| x.signum() * igamma(0.5, x*x) * sqrt_pi_inv).assert();
        
    should_the_same_mathfn(
        "erfc(x) = Γ(1/2, x^2)/√π",
            |x| erfc(x),
            |x| iGamma(0.5, x*x) * sqrt_pi_inv)
        .filter(|x| x >= 0.).assert();
}

//***** normal distribution *****
/// Return a value of the lower CDF (cumulative distribution function) of the normal distribution.
/// 
/// Ref: 『改訂新版 Cによる標準アルゴリズム事典』不完全ガンマ関数 (incomplete gamma function) igamma.c
pub fn p_normal(x: f64) -> f64 {
    if x >= 0. {
        0.5*(1. + p_gamma_normalizable(0.5, 0.5*x*x, LOG_PI_BY2))
    }else{
        0.5*q_gamma_normalizable(0.5, 0.5*x*x, LOG_PI_BY2)
    }
}

/// Return a value of the upper PDF of the normal distribution.
/// 
/// Ref: 『改訂新版 Cによる標準アルゴリズム事典』不完全ガンマ関数 (incomplete gamma function) igamma.c
pub fn q_normal(x: f64) -> f64 {
    if x >= 0. {
        0.5*q_gamma_normalizable(0.5, 0.5*x*x, LOG_PI_BY2)
    }else{
        0.5*(1. + p_gamma_normalizable(0.5, 0.5*x*x, LOG_PI_BY2))
    }
}

#[test]
fn test_the_normal_distribution_cdf_properties(){
    should_the_same_mathfn(
        "q_normal(x) = 1 - p_normal(x)", 
            |x| q_normal(x),
            |x| 1. - p_normal(x)).assert();

    use std::f64::consts::FRAC_1_SQRT_2;

    should_the_same_mathfn(
        "p_normal(x) = (1 + erf(x/√2))/2",
            |x| p_normal(x),
            |x| (1. + erf(x * FRAC_1_SQRT_2)) / 2.).assert();
}

//***** chi square distribution *****
/// Return a value of the lower CDF (cumulative distribution function) of the chi-square distribution.
/// 
/// Ref: 『改訂新版 Cによる標準アルゴリズム事典』不完全ガンマ関数 (incomplete gamma function) igamma.c
pub fn p_chi2(x: f64, n_freedom: f64) -> f64 {
    p_gamma(0.5 * n_freedom, 0.5 * x)
}

/// Return a value of the upper CDF of the chi-square distribution.
/// 
/// Ref: 『改訂新版 Cによる標準アルゴリズム事典』不完全ガンマ関数 (incomplete gamma function) igamma.c
pub fn q_chi2(x: f64, n_freedom: f64) -> f64 {
    q_gamma(0.5 * n_freedom, 0.5 * x)
}

#[test]
fn test_the_chi_square_distribution_cdf_properties(){
    should_the_same_mathfn2(
        "q_chi2(x, k) = 1 - p_chi2(x, k)",
            |x, k| q_chi2(x, k),
            |x, k| 1. - p_chi2(x, k))
        .var1(|v| v.name("k").range(1., 10.).is_integer(true).end()).assert();

    should_the_same_mathfn2(
        "p_chi2(x, k) = γ(k/2, x/2)/Γ(k/2)",
            |x, k| p_chi2(x, k),
            |x, k| p_gamma(k/2., x/2.))
        .var1(|v| v.name("k").range(1., 10.).is_integer(true).end()).assert();  
}

#[cfg(test)]
use super::test_util::*;
#[cfg(test)]
use super::gamma_fn::gamma;

#[test]
fn test_consts(){
    assert_eq!(LOG_PI_BY2, std::f64::consts::PI.ln() * 0.5);
}