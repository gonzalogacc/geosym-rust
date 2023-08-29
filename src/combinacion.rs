use crate::create_hs;

pub trait Model {
    fn generate_h(&self) -> Result<create_hs::SigmaH, create_hs::ModelError>;
}

#[allow(dead_code)]
#[derive(Debug, Copy, Clone)]
pub struct Powerlaw {
    pub m: usize,
    pub sigma: f32,
    pub kappa: f32,
    pub dt: f32, 
    pub units: create_hs::Units
}

#[allow(dead_code)]
#[derive(Debug, Copy, Clone)]
pub struct WhiteNoise {
    m: usize,
    sigma: f32,
}

impl Model for WhiteNoise {

    fn generate_h(&self) -> Result<create_hs::SigmaH, create_hs::ModelError> {
    // Create impulse function for White noise.
        create_hs::white(
            self.m, self.sigma
            )
    }
}

impl Model for Powerlaw {
    fn generate_h(&self) -> Result<create_hs::SigmaH, create_hs::ModelError> {
        create_hs::powerlaw(
            self.m, self.kappa, self.sigma, self.units, self.dt
            )
    }
}

// Gonzalo dice que un struct TotalNoise tiene que tener model_to_suma_de_todo
// como método
//
// Te dejo como tarea poner aca tu struct ideal.
//
// struct TotalNoise {
// }
//
pub fn model_to_h(ms: Vec<Box<dyn Model>> ) -> Vec<Vec<f32>> {
    ms.into_iter().map(|m| m.generate_h().unwrap().h).collect()
}

#[cfg(test)]
mod tests {
    use crate::combinacion::{
        Model, WhiteNoise, model_to_h
    };

    #[test]
    fn test_white() {
        let ms: Vec<Box<dyn Model>> = vec![ Box::new( WhiteNoise {m: 10, sigma: 1.0} ),
                       Box::new( WhiteNoise {m: 10, sigma: 1.0} ) ];
        assert_eq!( model_to_h(ms),vec![vec![1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],vec![1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]] );
    }
}
