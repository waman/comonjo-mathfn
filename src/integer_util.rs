pub trait IsInteger{
    fn is_integer(self) -> bool;
}

impl IsInteger for f64 {
    fn is_integer(self) -> bool {
        (self - self.round()).abs() < f64::EPSILON
    }
}