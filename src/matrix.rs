use core::{cmp::PartialEq, ops::{Mul, Add, Sub, Div, Deref, Neg}};
use std::{array::TryFromSliceError};
use std::mem::MaybeUninit;
use crate::permute::Permuter;
use num::traits::{Zero, One};

pub struct Matrix <T, const N: usize>{
    // TODO: store row or column order, convert only as needed
    m: [[T; N]; N]
}

pub trait Ring: Mul<Output = Self> + Div<Output = Self> + Sub<Output = Self> + Add<Output = Self> + Zero + One + Neg<Output = Self> {}

pub trait Determinant<T> {
    fn det(&self) -> T;
}

pub trait Inverse<T> : Determinant<T> {
    fn inv(&self) -> Self;
}

impl<T: Mul<Output = T> + Copy, const N: usize> Matrix<T, N> {

}

impl<T, const N: usize> Deref for Matrix<T, N> {
    type Target = [[T;N];N];

    fn deref(&self) -> &<Matrix<T,N> as Deref>::Target {
        &self.m
    }
}

impl<T: Copy, const N: usize> TryFrom<&[&[T]]> for Matrix<T, N> {
    type Error = TryFromSliceError;

    fn try_from(value: &[&[T]]) -> Result<Self, Self::Error> {
        unsafe {
            let mut tmp: MaybeUninit<[[T; N]; N]> = MaybeUninit::uninit();
            for (i, j) in value.iter().zip(0..N) {
                (*tmp.as_mut_ptr())[j] = <[T; N]>::try_from(*i)?;
            }
            Ok(Self{m: tmp.assume_init()})
        }
    }
}

impl<T: Copy, const N: usize> From<[[T; N]; N]> for Matrix<T, N> {
    fn from(value: [[T; N]; N]) -> Self {
        Self{m: value}
    }
}

impl<T: Mul<Output = T> + Add<Output = T> + Copy, const N: usize> Mul for &Matrix<T, N> {
    type Output = Matrix<T, N>;

    // TODO: sparse matrix multiplication?
    fn mul(self, rhs: Self) -> Self::Output {
        unsafe {
            let mut tmp: MaybeUninit<[[T; N]; N]> = MaybeUninit::uninit();
            for i in 0..N {
                for j in 0..N {
                    let mut sum = self[i][0] * rhs[0][j];
                    for k in 1..N {
                        sum = sum + self[i][k] * rhs[k][j]
                    }
                    (*tmp.as_mut_ptr())[i][j] = sum;
                }
            }
            Matrix{m: tmp.assume_init()}
        }
    }
}

impl<T: Mul<Output = T> + Copy, const N: usize> Mul<T> for Matrix<T, N> {
    type Output = Self;

    fn mul(self, rhs: T) -> Self::Output {
        let mut tmp = self.clone();
        for i in 0..N { for j in 0..N {
            tmp.m[i][j] = rhs * self.m[i][j];
        }}
        tmp
    }
}

impl<T: Div<Output = T> + Copy, const N: usize> Div<T> for Matrix<T, N> {
    type Output = Self;

    fn div(self, rhs: T) -> Self::Output {
        let mut tmp = self.clone();
        for i in 0..N { for j in 0..N {
            tmp.m[i][j] = self.m[i][j] / rhs;
        }}
        tmp
    }
}

impl<T: Copy, const N: usize> Clone for Matrix<T, N> {
    fn clone(&self) -> Self {
        Self{m: self.m.clone()}
    }
}

impl<T: PartialEq, const N: usize> PartialEq for Matrix<T, N> {

    fn eq(&self, other: &Self) -> bool {
        for i in 0..N {
            for j in 0..N {
                if self[i][j] != other[i][j] {return false};
            }
        }
        true
    }
}

impl<T: Ring + Copy, const N: usize> Determinant<T> for Matrix<T, N> {

    // TODO: fast Determinant?
    fn det(&self) -> T {
        match N {
            1 => self[0][0],
            2 => self[0][0] * self[1][1] - self[0][1] * self[1][0],
            _ => naive_det(&self.m),
        }
    }
}

fn naive_det<T, const N: usize>(m: &[[T; N]; N]) -> T
where T: Ring + Copy
{
    let mut permutations = Permuter::new(m.len());
    let mut sum = T::zero();
    let mut even_odd = T::one();
    while let Some(p) = permutations.next() {
        sum = sum + even_odd * p.iter().enumerate().map(|(i, j)| m[i][*j]).fold(T::zero(), |product, a| product * a);
        even_odd = -even_odd;
    }
    sum
}

impl<T: Ring + Copy, const N: usize> Inverse<T> for Matrix<T, N> {
    fn inv(&self) -> Self {
        match N {
            1 => (&[&[self[0][0] / self[0][0] / self[0][0]][0..]][0..]).try_into().expect(""),
            2 => TryInto::<Matrix<T, N>>::try_into(&[&[self[1][1], self[0][1].neg()][0..], &[self[1][0].neg(), self[0][0]][0..]][0..]).expect("") / self.det(),
            _ => unimplemented!(),
        }
    }
}

fn gray_code(n: u32) -> u32 {
    n ^ (n >> 1)
}

#[cfg(test)]
mod tests {
    use super::*;
    impl Ring for i32 {}

    #[test]
    fn test_eq() {
        let a: Matrix<i32, 2>= Into::<Matrix<i32, 2>>::into([[1, -1], [1, 1]]);
        let b: Matrix<i32, 2>= Into::<Matrix<i32, 2>>::into([[1, -1], [1, 1]]);
        assert!(a == b);
    }

    #[test]
    fn test_mul() {
        let a: Matrix<i32, 2>= Into::<Matrix<i32, 2>>::into([[1_i32, -1], [1, 1]]);
        let b: Matrix<i32, 2>= Into::<Matrix<i32, 2>>::into([[-1, -1], [2, 1]]);
        let c: Matrix<i32, 2>= Into::<Matrix<i32, 2>>::into([[-3, -2], [1, 0]]);
        assert!(&a * &b == c)
    }

    #[test]
    fn test_inv() {
        let a: Matrix<i32, 2>= Into::<Matrix<i32, 2>>::into([[2, 1], [1, 1]]);
        assert!(&a * &a.inv() == [[1,0],[0,1]].into());
        assert!(&a.inv() * &a == [[1,0],[0,1]].into());
        let a: Matrix<i32, 2>= Into::<Matrix<i32, 2>>::into([[-2, -7], [1, 4]]);
        assert!(&a * &a.inv() == [[1,0],[0,1]].into());
        assert!(&a.inv() * &a == [[1,0],[0,1]].into());
    }

    #[test]
    fn test_gray_code() {
        assert_eq!((0..4).map(|n| gray_code(n)).collect::<Vec<u32>>(), &[0,1,3,2])
    }
}
