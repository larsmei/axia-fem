use std::fmt;

#[derive(Debug)]
pub struct FemError(pub String);

impl fmt::Display for FemError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl std::error::Error for FemError {}

pub type Result<T> = std::result::Result<T, FemError>;

pub fn err<T>(msg: impl Into<String>) -> Result<T> {
    Err(FemError(msg.into()))
}
