// Esta funcion está super fuera de contexto pero la idea es que acá entre lo que haya que parsear
// y salga la unidad minima de especificación de modelo que Gonzalo quiera.
// Yo quería Box<dyn Model> pero el no. 

pub fn get_model(m: &str) -> Box<dyn Model> {
    match m {
        "white" => Box::new( WhiteNoise {m: 10, sigma: 1.0} ),
        "powerlaw" => Box::new( Powerlaw {m: 10, sigma: 1.0, kappa: 0.5, dt: 1.0, units: create_hs::Units::Mom } ),
        &_ => todo!()
    }
}

