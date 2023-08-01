
#[macro_export]
macro_rules! assert_eq_float {
    ($v1: expr, $v2: expr, $d: expr ) => {
        if ($v1 - $v2).abs() > $d {
            panic!("Values {:?} != {:?} differ more than {:?}", $v1, $v2, $d);
        }
    }
}

#[macro_export]
macro_rules! assert_eq_float_vec {
    ($v1: expr, $v2: expr, $d: expr ) => {
        if $v1.len() != $v2.len() {
            panic!("Vectors have different length {:?} : {:?}", $v1.len(), $v2.len());
        }
        for (i, value) in $v1.iter().enumerate(){
            println!("{:?}, {:?}", i, value);
            if (value - $v2[i]).abs() > $d {
                panic!("Values {:?} != {:?} differ more than {:?}", value, $v2[i], $d);
            }
        }
    }
}
