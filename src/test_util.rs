pub const EPS: f64 = 1e-12;
pub const BIG_VALUE: f64 = 1e13;

fn repeat<F>(n: usize, f: F)
    where F: Fn() -> ()
{
    for _ in 0..n { f(); }
}

pub fn assert_approximately(x: f64, y: f64, eps: f64, message: &str){
    let error = if y.abs() > 1.0 { y.abs() * eps }else{ eps };
    assert!((x - y).abs() < error, "{:?}: {} != {} +- {}", message, x, y, error);
}

pub fn assert_diverging(x: f64, big_value: f64, message: &str){
    assert!(x.is_nan() || x.is_infinite() || x.abs() >= big_value,
         "{:?}: {} is not diverging (big value: {})", message, x, big_value);
}

// pub struct FunctionTestProps<'a, F>
//     where F: Fn(f64) -> bool
// {
//     message: &'a str,
//     var_name: &'a str,
//     min: f64,
//     max: f64,
//     n: usize,
//     epsilon: f64,
//     filter: F
// }

// impl FunctionTestProps{
//     fn new<'a, F>(message: &str) -> FunctionTestProps<'a, F>
//         where F: Fn(f64) -> bool
//     {
//         FunctionTestProps{
//             message: message,
//             var_name: "x",
//             min: -5.0,
//             max: 5.0,
//             n: 1000,
//             epsilon: ESP,
//             filter: |x: f64| true
//         }
//     }
// }

// type FunctionTestProps = {
//     filter?: (x: number) => boolean;
//     message?: string;
//     varName?: string;
//     epsilon?: number;
//     n?: number;
//     min?: number;
//     max?: number;
// }

pub fn assert_the_same_mathfn<F, G, C>(f: F, g: G, filter: C, message: &str)
    where F: Fn(f64) -> f64, G: Fn(f64) -> f64, C: Fn(f64) -> bool
{
    let min = -5.0;
    let max: f64 = 5.0;
    let n: usize = 1000;
    let nf: f64 = n as f64;
    let delta = (max - min) / (n as f64);
    let error = EPS;

    fn assert_at<F, G, C>(x: f64, f: &F, g: &G, error: f64, filter: &C, message: &str)
        where F: Fn(f64) -> f64, G: Fn(f64) -> f64, C: Fn(f64) -> bool
    {
        if filter(x) {
            let yf = f(x);
            let yg = g(x);
            if both_are_nan(yf, yg) { return; }
            if both_are_infinite(yf, yg) { return; }

            if either_is_infinite(yf, yg) {
                assert!(yf.abs() > BIG_VALUE && yg.abs() > BIG_VALUE, "{}", message);
            }else{
                assert_approximately(f(x), g(x), error, 
                &format!("{} at x = {}", message, x));
            }
        }
    
        fn both_are_nan(x: f64, y: f64) -> bool {
            x.is_nan() && y.is_nan()
        }

        fn both_are_infinite(x: f64, y: f64) -> bool {
            x.is_infinite() && y.is_infinite()
        }

        fn either_is_infinite(x: f64, y: f64) -> bool {
            x.is_infinite() || y.is_infinite()
        }
    }

    let mut x = min;
    while x <= max {
        let x_round = (x * nf).round() / nf;
        assert_at(x_round, &f, &g, error, &filter, message);
        x += delta;
    }

    let interval = max - min;
    repeat(n, ||{
        let x = rand::random::<f64>() * interval + min;
        assert_at(x, &f, &g, error, &filter, message);
    });
}

// export function assertTheSameMathFunc(
//         f: (x: number) => number, g: (x: number) => number,
//         props: FunctionTestProps = {}){
//     const filter = props.filter ? props.filter : (x: number) => true;
//     const message = props.message ? props.message : '';
//     const varName = props.varName ? props.varName : 'x';
//     const epsilon = props.epsilon ? props.epsilon : EPS;
//     const n = props.n ? props.n : 1000;
//     const min = props.min ? props.min : -5;
//     const max = props.max ? props.max : 5;
// }