use core::{cmp::PartialEq, ops::{Mul, Add, Sub, Div, Deref, Neg}};
use std::{array::TryFromSliceError};
use std::mem::MaybeUninit;
use crate::permute::Permuter;
use num::traits::{Zero, One};

use std::rc::Rc;

#[derive(Clone)]
pub struct Matrix <T>{
    // TODO: store row or column order, convert only as needed
    m: Box<[T]>,
    dim: Dimension
}

#[derive(Debug, Clone, Copy)]
pub struct Dimension {
    row: usize,
    col: usize
}

pub struct SubMatrix<'a, T> {
    m: &'a [T],
    dim: ((usize, usize), Dimension)
}

impl<'a, T> SubMatrix<'a, T> {
    fn from_matrix<'b: 'a>(m: &'b Matrix<T>, dim: ((usize, usize), Dimension)) -> Self{
        SubMatrix{m: &m.m, dim}
    }
}

pub trait Ring: Mul<Output = Self> + Div<Output = Self> + Sub<Output = Self> + Add<Output = Self> + Zero + One + Neg<Output = Self> {}

impl<'a, 'b: 'a, T> From<&'b Matrix<T>> for SubMatrix<'a, T> {
    fn from(value: &'b Matrix<T>) -> Self {
        SubMatrix { m: &value.m, dim: ((0,0), value.dim) }
    }
}

impl<'a, T: Copy> SubMatrix<'a, T> {
    fn get(&self, loc: (usize, usize)) -> T{
        self.m[(self.dim.0.0 + loc.0) * self.dim.1.col + (self.dim.0.1 + loc.1)]
    }
}

impl<'a, T> SubMatrix<'a, T> {
    fn get_ref(&self, loc: (usize, usize)) -> &T{
        &self.m[(self.dim.0.0 + loc.0) * self.dim.1.col + (self.dim.0.1 + loc.1)]
    }
}

pub trait Determinant<T> {
    fn det(&self) -> T;
}

pub trait Inverse<T> : Determinant<T> {
    fn inv(&self) -> Matrix<T>;
}

impl<'a, T: Ring + Copy> SubMatrix<'a, T> {
    fn naive_det(&self) -> T {
        let mut permutations = Permuter::new(self.dim.1.row);
        let mut sum = T::zero();
        let mut even_odd = T::one();
        while let Some(p) = permutations.next() {
            sum = sum + even_odd * p.iter().enumerate().map(|(i, j)| self.get((i, *j))).fold(T::one(), |product, a| product * a);
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

// impl<T, const N: usize, const M: usize> Deref for Matrix<T, N, M> {
//     type Target = [[T;M];N];

//     fn deref(&self) -> &<Matrix<T,N,M> as Deref>::Target {
//         &self.m
//     }
// }

// impl<T: Copy, const N: usize, const M: usize> TryFrom<&[&[T]]> for Matrix<T, N, M> {
//     type Error = TryFromSliceError;

//     fn try_from(value: &[&[T]]) -> Result<Self, Self::Error> {
//         unsafe {
//             let mut tmp: MaybeUninit<[[T; M]; N]> = MaybeUninit::uninit();
//             for (i, j) in value.iter().zip(0..N) {
//                 (*tmp.as_mut_ptr())[j] = <[T; M]>::try_from(*i)?;
//             }
//             Ok(Self{m: tmp.assume_init()})
//         }
//     }
// }

impl<T: Copy, const N: usize, const M: usize> From<[[T; M]; N]> for Matrix<T> {
    fn from(value: [[T; M]; N]) -> Self {
        Self{m: value.into_iter().flatten().collect::<Vec<T>>().into_boxed_slice(), dim: Dimension{row: value.len(), col: value[0].len()}}
    }
}

impl<'a, 'b, T: Mul<Output = T> + Add<Output = T> + Copy> Mul<SubMatrix<'a, T>> for SubMatrix<'b, T>{
    type Output = Matrix<T>;

    // TODO: sparse matrix multiplication?
    fn mul(self, rhs: SubMatrix<'a, T>) -> Self::Output {
        assert!(self.dim.1.col == rhs.dim.1.row);
        let mut m = Vec::with_capacity(self.dim.1.row * rhs.dim.1.col);
        let pm: *mut T = m.as_mut_ptr();
        for i in 0..self.dim.1.row {
            for j in 0..rhs.dim.1.col {
                let mut sum = self.get((i, 0)) * rhs.get((0, j));
                for k in 1..self.dim.1.col {
                    sum = sum + self.get((i, k)) * rhs.get((k, j));
                }
                // This is safe becasue the access always occurs within m's allocated capacity
                unsafe {*(pm.add(i * rhs.dim.0.1 + j)) = sum; }
            }
        }
        Matrix{m: m.into_boxed_slice(), dim: Dimension{row: self.dim.1.row, col: rhs.dim.1.col}}
    }

}

impl<'a, 'b: 'a, 'c: 'a, T: Mul<Output = T> + Add<Output = T> + Copy> Mul<&'b Matrix<T>> for &'c Matrix<T>{
    type Output = Matrix<T>;

    fn mul(self, rhs: &'b Matrix<T>) -> Self::Output {
        Into::<SubMatrix<'a, T>>::into(self) * Into::<SubMatrix<'a, T>>::into(rhs)
    }
}

impl<'a, T: Mul<Output = T> + Copy> Mul<T> for SubMatrix<'a, T> {
    type Output = Matrix<T>;

    fn mul(self, rhs: T) -> Self::Output {
        let mut tmp = Vec::with_capacity(self.dim.1.row * self.dim.1.col);
        for i in 0..self.dim.1.row { for j in 0..self.dim.1.col {
            tmp.push(rhs * self.get((i, j)));
        }}
        Matrix{m: tmp.into_boxed_slice(), dim: self.dim.1}
    }
}

impl<'a, T: Div<Output = T> + Copy> Div<T> for SubMatrix<'a, T>{
    type Output = Matrix<T>;

    fn div(self, rhs: T) -> Self::Output {
        let mut tmp = Vec::with_capacity(self.dim.1.row * self.dim.1.col);
        for i in 0..self.dim.1.row { for j in 0..self.dim.1.col {
            tmp.push(self.get((i, j)) / rhs);
        }}
        Matrix{m: tmp.into_boxed_slice(), dim: self.dim.1}
    }
}

impl<T: Div<Output = T> + Copy> Div<T> for &Matrix<T> {
    type Output = Matrix<T>;
    fn div(self, rhs: T) -> Self::Output {
        Into::<SubMatrix<'_, T>>::into(self) / rhs
    }
}

impl<'a, T: PartialEq> PartialEq for SubMatrix<'a, T> {
    fn eq(&self, other: &Self) -> bool {
        assert!(self.dim.1.row == other.dim.1.row && self.dim.1.col == other.dim.1.col);
        for i in 0..self.dim.1.row {
            for j in 0..self.dim.1.col {
                if *self.get_ref((i, j)) != *other.get_ref((i, j)) {return false};
            }
        }
        true
    }
}

impl<T: PartialEq> PartialEq for Matrix<T> {
    fn eq(&self, other: &Self) -> bool {
        assert!(self.dim.row == other.dim.row && self.dim.col == other.dim.col);
        self.m.iter().zip(other.m.iter()).find(|(a, b)| a != b).is_none()
    }
}

impl<'a, T: Ring + Copy> Determinant<T> for SubMatrix<'a, T> {

    // TODO: fast Determinant?
    fn det(&self) -> T {
        assert!(self.dim.1.row == self.dim.1.col);
        match self.dim.1.row {
            1 => self.get((0, 0)),
            2 => self.get((0,0)) * self.get((1,1)) - self.get((0,1)) * self.get((1,0)),
            _ => self.naive_det(),
        }
    }
}

impl<T: Ring + Copy> Determinant<T> for Matrix<T> {
    fn det(&self) -> T {
        Into::<SubMatrix<'_, T>>::into(self).det()
    }
}

impl<'a, T: Ring + Copy> Inverse<T> for SubMatrix<'a, T> {

    fn inv(&self) -> Matrix<T> {
        match self.dim.1.row {
            1 => ([[T::one() / self.get((0,0))]]).try_into().expect(""),
            2 => &Into::<Matrix<T>>::into([[self.get((1,1)), -self.get((0,1))], [-self.get((1,0)), self.get((0,0))]]) / self.det(),
            _ => unimplemented!(),
        }
    }
}

impl<T: Ring + Copy> Inverse<T> for Matrix<T> {
    fn inv(&self) -> Matrix<T> {
        Into::<SubMatrix<'_, T>>::into(self).inv()
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
        let a: Matrix<i32> = [[1, -1], [1, 1]].into();
        let b = [[1, -1], [1, 1]].into();
        assert!(a == b);
    }

    #[test]
    fn test_mul() {
        let a: Matrix<i32> = [[1_i32, -1], [1, 1]].into();
        let b = [[-1, -1], [2, 1]].into();
        let c = [[-3, -2], [1, 0]].into();
        let d = [[3,4,5], [-1,0,1]].into();
        let e = [[4, 4, 4], [2, 4, 6]].into();
        assert!(&a * &b == c);
        assert!(&a * &d == e);
    }


    #[test]
    fn test_det() {
        let a: Matrix<i32> = [[1,2,3],[-1,-2,-3],[1,0,-1]].into();
        let b: Matrix<i32> = [[1,2,3],[-3,-2,-1],[1,0,-1]].into();
        let c: Matrix<i32> = [[1,2,3],[-3,-2,-1],[1,0,1]].into();
        assert_eq!(a.det(), 0);
        assert_eq!(b.det(), 0);
        assert_eq!(c.det(), 8);
    }

    #[test]
    fn test_inv_base() {
        let a: Matrix<i32> = [[2, 1], [1, 1]].into();
        assert!(&a * &a.inv() == [[1,0],[0,1]].into());
        assert!(&a.inv() * &a == [[1,0],[0,1]].into());
        let a: Matrix<i32> = [[-2, -7], [1, 4]].into();
        assert!(&a * &a.inv() == [[1,0],[0,1]].into());
        assert!(&a.inv() * &a == [[1,0],[0,1]].into());
    }

    // #[test]
    // fn test_gray_code() {
    //     assert_eq!((0..4).map(|n| gray_code(n)).collect::<Vec<u32>>(), &[0,1,3,2])
    // }
}
