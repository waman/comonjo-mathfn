use comonjo_mathfn::gamma;

fn main(){
    for x in -3..=5 {
        println!("Γ({}) = {}", x, gamma(x as f64));
    }
    println!("Γ(1/2) = {} (~ √π = {})", gamma(0.5), std::f64::consts::PI.sqrt());
}