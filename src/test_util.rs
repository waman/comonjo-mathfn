use rand::random;

pub const EPS: f64 = 1e-12;
pub const BIG_VALUE: f64 = 1e13;

pub trait IsInteger {
    fn is_integer(&self, delta: f64) -> bool;
    fn is_non_positive_integer(&self, delta: f64) -> bool;
}

impl IsInteger for f64 {
    fn is_integer(&self, delta: f64) -> bool {
        (*self - self.round()).abs() <= delta
    }
    
    fn is_non_positive_integer(&self, delta: f64) -> bool {
        *self <= delta && self.is_integer(delta)
    }
}

pub fn repeat<F>(n: usize, f: F)
    where F: Fn() -> ()
{
    for _ in 0..n { f(); }
}

pub fn rand(min: f64, max: f64) -> f64 {
    min + (max - min) * random::<f64>()
}

pub fn assert_approximately(x: f64, y: f64, eps: f64, message: &str){
    if x.is_nan() || y.is_nan() {
         assert!(x.is_nan() && y.is_nan(),
            "{} --- only either value is nan: {}, {}", message, x, y);
    
    }else{
        let x_abs = x.abs();
        let y_abs = y.abs();
        
        if x_abs >= BIG_VALUE || y_abs >= BIG_VALUE {
            assert!(x * y > 0. &&  x_abs >= BIG_VALUE && y_abs >= BIG_VALUE,
                "{} --- only either value diverge: {}, {}", message, x, y);

        }else{
            if y_abs > 1. {
                let error = y_abs * eps;
                assert_approx(x, y, error, message);
            } else {
                assert_approx(x, y, eps, message);
            }
        }
    }

    fn assert_approx(x: f64, y: f64, error: f64, message: &str){
        assert!((x - y).abs() < error, "{}", format!("{}: {} != {} +- {} (error ratio: 1e-{})",
                message, x, y, error, ratio(x, y)));
    }

    fn ratio(x: f64, y: f64) -> usize { 
        let r = ((x - y) / x).abs();
        let mut s = 0.1;
        let mut n = 1;
        while s > r {
            s /= 10.;
            n += 1;
        }
        n
    } 
}

pub fn assert_diverging(x: f64, big_value: f64, message: &str){
    assert!(x.is_nan() || x.is_infinite() || x.abs() >= big_value,
         "{:?}: {} is not diverging (big value: {})", message, x, big_value);
}

pub struct TestVariable<'a>{
    _name: &'a str,
    _range: (f64, f64),
    _n: usize,
    _is_integer: bool
}

impl<'a> TestVariable<'a> {

    pub fn name(&mut self, new_name: &'a str) -> &mut Self {
        self._name = new_name;
        self
    }

    pub fn range(&mut self, min: f64, max: f64) -> &mut Self {
        self._range = (min, max);
        self
    }

    pub fn is_integer(&mut self, b: bool) -> &mut Self {
        self._is_integer = b;
        self
    }

    pub fn end(&self){()}
    
    fn at_even_intervals<F>(&self, f: F) where F: Fn(f64) {
        if self._is_integer {
            let min = self._range.0.round() as usize;
            let max = self._range.1.round() as usize;
            for i in min..=max {
                f(i as f64);
            }

        } else {
            let (min, max) = self._range; 
            let nf = self._n as f64;
            let delta = (max - min) / nf;

            let mut x = min;
            while x <= max {
                let xr = (x * nf).round() / nf;
                f(xr); 
                x += delta;
            }
        }
    }

    fn at_random<F>(&self, f: F)  where F: Fn(f64) {
        if self._is_integer {
            self.at_even_intervals(f);
        } else {
            let (min, max) = self._range;
            let interval = max - min;
            repeat(self._n, || {
                let x = min + interval * rand::random::<f64>();
                f(x);
            });
        }
    }
}

pub struct MathFnTest<'a, F, G>
    where F: Fn(f64) -> f64 + 'a,
          G: Fn(f64) -> f64 + 'a
{
    message: &'a str,
    f: F,
    g: G,
    _var: TestVariable<'a>,
    _filter: Box<dyn Fn(f64) -> bool>,
    _epsilon: f64
}

impl<'a, F, G> MathFnTest<'a, F, G>
    where F: Fn(f64) -> f64 + 'a,
          G: Fn(f64) -> f64 + 'a
{
    fn new<'t>(message: &'t str, f: F, g: G) -> MathFnTest<'t, F, G>{
        MathFnTest{
            message, f, g,
            _var: TestVariable{ _name: "x", _range: (-5., 5.), _n: 1000, _is_integer: false },
            _filter: Box::new(|_| true),
            _epsilon: EPS
        }
    }

    pub fn filter<P>(&mut self, pred: P) -> &mut Self
        where P: Fn(f64) -> bool + 'static
    {
        self._filter = Box::new(pred);
        self
    }

    pub fn var0<S>(&mut self, settings: S) -> &mut Self
        where S: Fn(&mut TestVariable<'a>)
    {
        settings(&mut self._var);
        self
    }

    pub fn epsilon(&mut self, eps: f64) -> &mut Self {
        self._epsilon = eps;
        self
    }

    fn assert_at(&self, x: f64){
        if (self._filter)(x) {
            assert_approximately((self.f)(x), (self.g)(x), self._epsilon,
                &format!("{} at {} = {}", self.message, self._var._name, x));
        }
    }
    
    pub fn assert(&self){
        self._var.at_even_intervals(|x| {
            self.assert_at(x);
        });
        
        if self._var._is_integer { return; }
        self._var.at_random(|x| {
            self.assert_at(x);
        });
    }
}

pub fn should_the_same_mathfn<'a, F, G>(message: &'a str, f: F, g: G) -> MathFnTest<'a, F, G> 
    where F: Fn(f64) -> f64 + 'a,
          G: Fn(f64) -> f64 + 'a
{
    MathFnTest::new(message, f, g)
}

pub struct MathFnTest2<'a, F, G>
    where F: Fn(f64, f64) -> f64 + 'a,
          G: Fn(f64, f64) -> f64 + 'a
{
    message: &'a str,
    f: F,
    g: G,
    _var0: TestVariable<'a>,
    _var1: TestVariable<'a>,
    _filter: Box<dyn Fn(f64, f64) -> bool>,
    _epsilon: f64
}

impl<'a, F, G> MathFnTest2<'a, F, G>
    where F: Fn(f64, f64) -> f64 + 'a,
          G: Fn(f64, f64) -> f64 + 'a
{
    fn new<'t>(message: &'t str, f: F, g: G) -> MathFnTest2<'t, F, G>{
        let v = TestVariable{ _name: "y", _range: (-5., 5.), _n: 100, _is_integer: false };
        MathFnTest2{
            message, f, g,
            _var0: TestVariable{ _name: "x", .. v },
            _var1: v,
            _filter: Box::new(|_x, _y| true),
            _epsilon: EPS
        }
    }

    pub fn filter<P>(&mut self, pred: P) -> &mut Self
        where P: Fn(f64, f64) -> bool + 'static
    {
        self._filter = Box::new(pred);
        self
    }

    pub fn var0<S>(&mut self, settings: S) -> &mut Self
        where S: Fn(&mut TestVariable<'a>)
    {
        settings(&mut self._var0);
        self
    }

    pub fn var1<'b, S>(&mut self, settings: S) -> &mut Self
        where S: Fn(&mut TestVariable<'a>)
    {
        settings(&mut self._var1);
        self
    }

    pub fn epsilon(&mut self, eps: f64) -> &mut Self {
        self._epsilon = eps;
        self
    }

    fn assert_at(&self, x: f64, y: f64){
        if (self._filter)(x, y) {
            assert_approximately((self.f)(x, y), (self.g)(x, y), self._epsilon,
                &format!("{} at ({}, {}) = ({}, {})",
                    self.message, self._var0._name, self._var1._name, x, y));
        }
    }
    
    pub fn assert(&self){
        self._var1.at_even_intervals(|y| {
            self._var0.at_even_intervals(|x| {
                self.assert_at(x, y);
            });
        });
        
        if self._var0._is_integer && self._var1._is_integer { return; }
        self._var1.at_random(|y| {
            self._var0.at_random(|x| {
                self.assert_at(x, y);
            });
        });
    }
}

pub fn should_the_same_mathfn2<'a, F, G>(message: &'a str, f: F, g: G) -> MathFnTest2<'a, F, G> 
    where F: Fn(f64, f64) -> f64 + 'static,
          G: Fn(f64, f64) -> f64 + 'static
{
    MathFnTest2::new(message, f, g)
}