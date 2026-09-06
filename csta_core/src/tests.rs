use crate::{geometry::*, vec2::*, vec3::*, vec4::*};
macro_rules! vector_test {
    ($name:ident,$a:ident,$b:ident,$n:expr) => {
        #[test]
        fn $name() {
            let values: [f64; $n] = std::array::from_fn(|i| i as f64 + 1.0);
            let a = $a::from(values);
            let b = $b::from(a);
            let round = $a::from(b);
            assert_eq!(a, round);
            assert_eq!(<[f64; $n]>::from(&a), values);
            assert_eq!(&a + &a, a * 2.0);
            assert_eq!(&a - &a, $a::default());
            assert_eq!(a.dot(&a), values.iter().map(|x| x * x).sum::<f64>());
            assert_eq!(a.distance(&a), 0.0);
            assert!((a.normalize().len() - 1.0).abs() < 1e-14);
            for divisor in [0.0, -0.0, 1e-310, f64::INFINITY, f64::NAN, 2.0] {
                let x: [f64; $n] = (a / divisor).into();
                let mut y = a;
                y /= divisor;
                let y: [f64; $n] = y.into();
                for (x, y) in x.into_iter().zip(y) {
                    assert!(x == y || (x.is_nan() && y.is_nan()));
                }
            }
            for scale in [1e-310, 1e308] {
                assert!((($a::from([scale; $n])).normalize().len() - 1.0).abs() < 1e-14);
            }
            for scale in [1e-40_f32, 1e38] {
                assert!((($b::from([scale; $n])).normalize().len() - 1.0).abs() < 1e-6);
            }
            assert!($a::default().try_normalize().is_none());
            assert!($a::from([f64::NAN; $n]).try_normalize().is_none());
            assert!($b::from([f32::INFINITY; $n]).try_normalize().is_none());
        }
    };
}
vector_test!(vec2, Vec2f64, Vec2f32, 2);
vector_test!(vec3, Vec3f64, Vec3f32, 3);
vector_test!(vec4, Vec4f64, Vec4f32, 4);
#[test]
fn geometry() {
    let b = PeriodicBox::new([10.0, 20.0]).unwrap();
    assert_eq!(b.wrap([-1.0, 61.0]).unwrap(), [9.0, 1.0]);
    assert_eq!(
        b.displacement([9.0, 0.0], [1.0, 10.0]).unwrap(),
        [2.0, -10.0]
    );
    for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(PeriodicBox::new([bad]).is_err());
    }
    assert!(b.wrap([f64::NAN, 0.0]).is_err());
    for boundary in [Boundary::Open, Boundary::Periodic] {
        let l = Lattice2D::new(3, 4, boundary).unwrap();
        assert_eq!(
            l.bonds().count(),
            if boundary == Boundary::Open { 17 } else { 24 }
        );
        for (i, j) in l.bonds() {
            assert!(l.neighbors(j).unwrap().contains(&i));
        }
    }
    assert!(Lattice2D::new(2, 3, Boundary::Periodic).is_err());
    assert!(Lattice2D::new(usize::MAX, 2, Boundary::Open).is_err());
    assert_eq!(
        Lattice2D::new(1, 1, Boundary::Open)
            .unwrap()
            .bonds()
            .count(),
        0
    );
}

#[cfg(feature = "serde")]
#[test]
fn serde_vector_roundtrips() {
    macro_rules! roundtrip {
        ($t:ty,$value:expr) => {{
            let value: $t = $value;
            let json = serde_json::to_vec(&value).unwrap();
            let result: $t = serde_json::from_slice(&json).unwrap();
            assert_eq!(value, result);
        }};
    }
    roundtrip!(Vec2f32, Vec2f32(1.0, 2.0));
    roundtrip!(Vec2f64, Vec2f64(1.0, 2.0));
    roundtrip!(Vec3f32, Vec3f32(1.0, 2.0, 3.0));
    roundtrip!(Vec3f64, Vec3f64(1.0, 2.0, 3.0));
    roundtrip!(Vec4f32, Vec4f32(1.0, 2.0, 3.0, 4.0));
    roundtrip!(Vec4f64, Vec4f64(1.0, 2.0, 3.0, 4.0));
}
#[test]
fn periodic_translation_and_brute_force_images() {
    let cell = PeriodicBox::new([3.0, 5.0]).unwrap();
    for x in [-100.1, -3.0, -0.1, 0.0, 0.1, 2.9, 30.1] {
        let p = cell.wrap([x, x]).unwrap();
        let shifted = cell.wrap([x + 300.0, x - 500.0]).unwrap();
        for i in 0..2 {
            assert!((p[i] - shifted[i]).abs() < 1e-12);
        }
        let r = cell.displacement([x, 0.2], [1.1, 3.2]).unwrap();
        let a = cell.wrap([x, 0.2]).unwrap();
        let mut best = f64::INFINITY;
        for i in -1..=1 {
            for j in -1..=1 {
                let d = [1.1 - a[0] + 3.0 * i as f64, 3.2 - a[1] + 5.0 * j as f64];
                best = best.min(d[0] * d[0] + d[1] * d[1]);
            }
        }
        assert!((r[0] * r[0] + r[1] * r[1] - best).abs() < 1e-12);
    }
}
#[test]
fn tuples_accessors_and_f32_operators() {
    assert_eq!(<(f64, f64)>::from(Vec2f64::new(1.0, 2.0)), (1.0, 2.0));
    assert_eq!(
        <(f64, f64, f64)>::from(Vec3f64::new(1.0, 2.0, 3.0)),
        (1.0, 2.0, 3.0)
    );
    assert_eq!(
        <(f64, f64, f64, f64)>::from(Vec4f64::new(1.0, 2.0, 3.0, 4.0)),
        (1.0, 2.0, 3.0, 4.0)
    );
    let v = Vec4f32::from((1.0, 2.0, 3.0, 4.0));
    assert_eq!([v.x(), v.y(), v.z(), v.w()], [1.0, 2.0, 3.0, 4.0]);
    macro_rules! check {
        ($t:ident,$n:expr) => {{
            let v = $t::from([2.0_f32; $n]);
            let mut a = v;
            a += &v;
            assert_eq!(a, v * 2.0);
            a -= &v;
            assert_eq!(a, v);
            a *= 2.0;
            assert_eq!(a, 2.0 * v);
            assert_eq!(-a, $t::from([-4.0; $n]));
            for s in [0.0, -0.0, 1e-40, f32::INFINITY, f32::NAN] {
                let x: [f32; $n] = (v / s).into();
                let mut y = v;
                y /= s;
                let y: [f32; $n] = y.into();
                for (x, y) in x.into_iter().zip(y) {
                    assert!(x == y || (x.is_nan() && y.is_nan()));
                }
            }
        }};
    }
    check!(Vec2f32, 2);
    check!(Vec3f32, 3);
    check!(Vec4f32, 4);
}
#[cfg(feature = "serde")]
#[test]
fn deserialized_geometry_validates_constructors() {
    assert!(serde_json::from_str::<PeriodicBox<1>>(r#"{"lengths":[0.0]}"#).is_err());
    assert!(
        serde_json::from_str::<Lattice2D>(r#"{"width":0,"height":3,"boundary":"Open"}"#).is_err()
    );
    let b = PeriodicBox::new([2.0; 3]).unwrap();
    let decoded: PeriodicBox<3> = serde_json::from_slice(&serde_json::to_vec(&b).unwrap()).unwrap();
    assert_eq!(b, decoded);
}
#[cfg(feature = "checkpoint")]
#[test]
fn relative_checkpoint_path() {
    let file = format!(".csta-relative-test-{}.json", std::process::id());
    crate::checkpoint::save(&file, "scalar", &1.0_f64).unwrap();
    assert_eq!(
        crate::checkpoint::load::<f64>(&file, "scalar").unwrap(),
        1.0
    );
    std::fs::remove_file(file).unwrap();
}
