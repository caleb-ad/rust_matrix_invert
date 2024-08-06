use core::{cmp::PartialEq, ops::{Mul, Add, Sub, Div, Deref, Neg}};
use std::{array::TryFromSliceError};
use std::mem::MaybeUninit;
use crate::permute::Permuter;
use num::traits::{Zero, One};

#[derive(Clone, Copy)]
pub struct Matrix <T, const N: usize, const M: usize>{
    // TODO: store row or column order, convert only as needed
    m: [[T; M]; N]
}

pub trait Ring: Mul<Output = Self> + Div<Output = Self> + Sub<Output = Self> + Add<Output = Self> + Zero + One + Neg<Output = Self> {}

pub trait Determinant<T> {
    fn det(&self) -> T;
}

pub trait Inverse<T> : Determinant<T> {
    fn inv(&self) -> Self;
}

impl<T: Ring + Copy, const N: usize> Matrix<T, N, N> {
    fn naive_det(&self) -> T
    where T: Ring + Copy
    {
        let mut permutations = Permuter::new(self.m.len());
        let mut sum = T::zero();
        let mut even_odd = T::one();
        while let Some(p) = permutations.next() {
            sum = sum + even_odd * p.iter().enumerate().map(|(i, j)| self.m[i][*j]).fold(T::zero(), |product, a| product * a);
            even_odd = -even_odd;
        }
        sum
    }

    fn invert_by_block(&self) -> Self{
        match N {
            2 | 3 => self.inv(),
            _ => {
                let A: Matrix<T, {N/2}, N/2> = self.m[0..N/2].iter().map(|a| &a[0..N/2]).collect().into();
                let B: Matrix<T, N/2, (N+1)/2> = self.m[0..N/2].iter().map(|a| &a[N/2..]).collect().into();
                let D: Matrix<T, (N + 1)/2, (N + 1)/2> = self.m.iter().map(|a| &a[N/2..]).collect().into();


                asdf
            }
        }
    }

}

impl<T, const N: usize, const M: usize> Deref for Matrix<T, N, M> {
    type Target = [[T;M];N];

    fn deref(&self) -> &<Matrix<T,N,M> as Deref>::Target {
        &self.m
    }
}

impl<T: Copy, const N: usize, const M: usize> TryFrom<&[&[T]]> for Matrix<T, N, M> {
    type Error = TryFromSliceError;

    fn try_from(value: &[&[T]]) -> Result<Self, Self::Error> {
        unsafe {
            let mut tmp: MaybeUninit<[[T; M]; N]> = MaybeUninit::uninit();
            for (i, j) in value.iter().zip(0..N) {
                (*tmp.as_mut_ptr())[j] = <[T; M]>::try_from(*i)?;
            }
            Ok(Self{m: tmp.assume_init()})
        }
    }
}

impl<T: Copy, const N: usize, const M: usize> From<[[T; M]; N]> for Matrix<T, N, M> {
    fn from(value: [[T; M]; N]) -> Self {
        Self{m: value}
    }
}

impl<T: Mul<Output = T> + Add<Output = T> + Copy, const N: usize, const M: usize, const L: usize> Mul<&Matrix<T,M,L>> for &Matrix<T, N, M> {
    type Output = Matrix<T, N, L>;

    // TODO: sparse matrix multiplication?
    fn mul(self, rhs: &Matrix<T,M,L>) -> Self::Output {
        unsafe {
            let mut tmp: MaybeUninit<[[T; L]; N]> = MaybeUninit::uninit();
            for i in 0..N {
                for j in 0..L {
                    let mut sum = self[i][0] * rhs[0][j];
                    for k in 1..M {
                        sum = sum + self[i][k] * rhs[k][j]
                    }
                    (*tmp.as_mut_ptr())[i][j] = sum;
                }
            }
            Matrix{m: tmp.assume_init()}
        }
    }
}

impl<T: Mul<Output = T> + Copy, const N: usize, const M: usize> Mul<T> for Matrix<T, N, M> {
    type Output = Self;

    fn mul(self, rhs: T) -> Self::Output {
        let mut tmp = self.clone();
        for i in 0..N { for j in 0..N {
            tmp.m[i][j] = rhs * self.m[i][j];
        }}
        tmp.into()
    }
}

impl<T: Div<Output = T> + Copy, const N: usize, const M: usize> Div<T> for Matrix<T, N, M> {
    type Output = Self;

    fn div(self, rhs: T) -> Self::Output {
        let mut tmp = self.clone();
        for i in 0..N { for j in 0..N {
            tmp.m[i][j] = self.m[i][j] / rhs;
        }}
        tmp.into()
    }
}

impl<T: PartialEq, const N: usize, const M: usize> PartialEq for Matrix<T, N, M> {

    fn eq(&self, other: &Self) -> bool {
        for i in 0..N {
            for j in 0..N {
                if self[i][j] != other[i][j] {return false};
            }
        }
        true
    }
}

impl<T: Ring + Copy, const N: usize> Determinant<T> for Matrix<T, N, N> {

    // TODO: fast Determinant?
    fn det(&self) -> T {
        match N {
            1 => self[0][0],
            2 => self[0][0] * self[1][1] - self[0][1] * self[1][0],
            _ => self.naive_det(),
        }
    }
}

impl<T: Ring + Copy, const N: usize> Inverse<T> for Matrix<T, N, N> {
    fn inv(&self) -> Self {
        match N {
            1 => (&[&[self[0][0] / self[0][0] / self[0][0]][0..]][0..]).try_into().expect(""),
            2 => TryInto::<Matrix<T, N, N>>::try_into(&[&[self[1][1], self[0][1].neg()][0..], &[self[1][0].neg(), self[0][0]][0..]][0..]).expect("") / self.det(),
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
        let a = Into::<Matrix<i32, 2,2>>::into([[1, -1], [1, 1]]);
        let b = Into::<Matrix<i32, 2,2>>::into([[1, -1], [1, 1]]);
        assert!(a == b);
    }

    #[test]
    fn test_mul() {
        let a = Into::<Matrix<i32, 2,2>>::into([[1_i32, -1], [1, 1]]);
        let b = Into::<Matrix<i32, 2,2>>::into([[-1, -1], [2, 1]]);
        let c = Into::<Matrix<i32, 2,2>>::into([[-3, -2], [1, 0]]);
        let d = Into::<Matrix<i32, 2, 3>>::into([[3,4,5], [-1,0,1]]);
        let e = Into::<Matrix<i32, 2, 3>>::into([[4, 4, 4], [2, 4, 6]]);
        assert!(&a * &b == c);
        assert!(&a * &d == e);
    }

    #[test]
    fn test_inv() {
        let a = Into::<Matrix<i32, 2,2>>::into([[2, 1], [1, 1]]);
        assert!(&a * &a.inv() == [[1,0],[0,1]].into());
        assert!(&a.inv() * &a == [[1,0],[0,1]].into());
        let a = Into::<Matrix<i32, 2,2>>::into([[-2, -7], [1, 4]]);
        assert!(&a * &a.inv() == [[1,0],[0,1]].into());
        assert!(&a.inv() * &a == [[1,0],[0,1]].into());
    }

    #[test]
    fn test_gray_code() {
        assert_eq!((0..4).map(|n| gray_code(n)).collect::<Vec<u32>>(), &[0,1,3,2])
    }
}
