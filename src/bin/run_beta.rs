use comonjo_mathfn::beta;

fn main(){
    for i in 0..=6 {
        for j in 0..=6 {
            let x = (i as f64) - 3.;
            let y = (j as f64) - 3.;
            let tab = " ".repeat(j);
            println!("{}Β({}, {}) = {},  ", tab, x, y, beta(x, y));
        }
    }
    println!("Β(1/2, 1/2) = {} (~ π)", beta(0.5, 0.5));
}