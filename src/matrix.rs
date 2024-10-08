use core::{cmp::PartialEq, ops::{Mul, Add, Sub, Div, Deref, Neg}};
use std::{marker::Copy};
use crate::{graph::{Graph, Undirected, Cycle}, permute::Permuter};
use num::{complex::ComplexFloat, traits::{One, Zero}};
use std::str::FromStr;

#[derive(Debug, Clone)]
pub struct Matrix <T>{
    // TODO: store row or column order, convert only as needed
    m: Box<[T]>,
    pub dim: Dimension
}

impl<T: Copy> Matrix<T> {

    // The caller is responsible for initializing every T
    unsafe fn new_uninit(dim: Dimension) -> Self {
        Matrix {
            m: Box::from_raw(
                std::slice::from_raw_parts_mut(
                    std::alloc::alloc(std::alloc::Layout::array::<T>(dim.col * dim.row).expect("matrix too large")) as *mut T
                    , dim.col * dim.row) as *mut [T]
            ),
            dim
        }
    }

    fn from_submatrices(a: SubMatrix<'_, T>, b: SubMatrix<'_, T>, c: SubMatrix<'_, T>, d: SubMatrix<'_, T>) -> Matrix<T> {
        assert!(c.dim.1.row == d.dim.1.row && b.dim.1.col == d.dim.1.col && a.dim.1.col == c.dim.1.col && a.dim.1.row == b.dim.1.row);
        let dim = Dimension{row: a.dim.1.col + b.dim.1.col, col: a.dim.1.row + c.dim.1.row};
        unsafe {
            let mut tmp = Matrix::new_uninit(dim);
            tmp.copy_submatrix(b, (0, a.dim.1.col));
            tmp.copy_submatrix(c, (a.dim.1.row, 0));
            tmp.copy_submatrix(d, (a.dim.1.row, a.dim.1.col));
            tmp.copy_submatrix(a, (0,0));
            tmp
        }
    }

    // copy submatrix into self, submatrix should not be larger than self
    fn copy_submatrix(&mut self, a: SubMatrix<'_, T>, base: (usize, usize)) {
        for i in 0..a.dim.1.row {
            for j in 0..a.dim.1.col {
                self.m[(base.0 + i) * self.dim.col + base.1 + j] = a.get((i,j))
            }
        }
    }

}

impl<T: Clone + One + Zero> Matrix<T> {
    fn zero(n: usize, m: usize) -> Self {
        let mut tmp = Vec::with_capacity(n*m);
        tmp.resize(n*m, T::zero());
        Matrix{m: tmp.into_boxed_slice(), dim: Dimension{row: n, col: m}}
    }

    fn identity(n: usize) -> Self {
        let mut tmp = Matrix::zero(n, n);
        for i in 0..n {tmp.m[i*n + i] = T::one()}
        tmp
    }
}

impl<T> Matrix<T> {
    pub fn get_ref(&self, loc: (usize, usize)) -> &T {
        &self.m[loc.0 * self.dim.col + loc.1]
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Dimension {
    pub row: usize,
    pub col: usize
}

#[derive(Clone, Copy)]
pub struct SubMatrix<'a, T> {
    m: &'a Matrix<T>,
    pub dim: ((usize, usize), Dimension)
}

impl<'a, T> SubMatrix<'a, T> {
    fn from_matrix<'b: 'a>(m: &'b Matrix<T>, dim: ((usize, usize), Dimension)) -> Self {
        SubMatrix{m, dim}
    }

    fn from_submatrix(m: &SubMatrix<'a, T>, dim: ((usize, usize), Dimension)) -> Self {
        SubMatrix{m: m.m, dim: ((m.dim.0.0 + dim.0.0, m.dim.0.1 + dim.0.1), dim.1)}
    }

    pub fn get_ref(&self, loc: (usize, usize)) -> &T{
        &self.m.m[(self.dim.0.0 + loc.0) * self.m.dim.col + (self.dim.0.1 + loc.1)]
    }
}

impl<'a, T: Copy> SubMatrix<'a, T> {
    fn get(&self, loc: (usize, usize)) -> T{
        self.m.m[(self.dim.0.0 + loc.0) * self.m.dim.col + (self.dim.0.1 + loc.1)]
    }
}

impl<'a> SubMatrix<'a, f64> {
    fn nearly_equal(&self, other:  &SubMatrix<'_, f64>) -> bool {
        const THRESHHOLD: f64 = 0.000001;
        assert!(self.dim == other.dim);
        for i in 0..self.dim.1.row {
            for j in 0..self.dim.1.col {
                if self.get((i,j)) - other.get((i,j)) > THRESHHOLD {return false;}
            }
        }
        true
    }
}

impl<'a, T: Ring + Copy + PartialEq> SubMatrix<'a, T> {
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

    fn invert_3x3(&self) -> Matrix<T> {
        let cross = |i,j,k,l| self.get((i,j))*self.get((k,l))-self.get((i,l))*self.get((k,j));
        let tmp = [
            [cross(1,1,2,2), -cross(0,1,2,2), cross(0,1,1,2)],
            [-cross(1,0,2,2), cross(0,0,2,2), -cross(0,0,1,2)],
            [cross(1,0, 2, 1), -cross(0,0,2,1),cross(0,0,1,1)]
        ];
        Into::<Matrix<T>>::into(tmp) / self.det()
    }

    fn invert_by_block(self) -> Matrix<T> {
        let N = self.dim.1.row;
        match N {
            1 | 2 | 3 => self.inv(),
            _ => {
                let A = SubMatrix::from_submatrix(&self, ((0,0), Dimension{row: N/2, col: N/2}));
                let B = SubMatrix::from_submatrix(&self, ((0, N/2), Dimension{row: N/2, col: (N+1)/2}));
                let C = SubMatrix::from_submatrix(&self, ((N/2, 0), Dimension{row: (N+1)/2, col: N/2}));
                let D = SubMatrix::from_submatrix(&self, ((N/2, N/2), Dimension{row: (N+1)/2, col: (N+1)/2}));

                let t1: Matrix<T> = SubMatrix::invert_by_block((&(A - (B * D.invert_by_block() * C))).into());
                let t2: Matrix<T> = SubMatrix::invert_by_block((&(D - (C * A.invert_by_block() * B))).into());
                let t3: Matrix<T> =  (B * D.invert_by_block()) * -T::one();
                let t4: Matrix<T> =  (C * A.invert_by_block()) * -T::one();

                Matrix::from_submatrices(
                    (&t1).into(),
                    (&Matrix::zero(N/2,(N+1)/2)).into(),
                    (&Matrix::zero((N+1)/2, N/2)).into(),
                    (&t2).into())
                * Matrix::from_submatrices(
                    (&Matrix::identity(N/2)).into(),
                    (&t3).into(),
                    (&t4).into(),
                    (&Matrix::identity((N+1)/2)).into())
            }
        }
    }

    fn invert_sparse(self) -> Matrix<T> {
        const SPARSITY_THRESHHOLD: usize = 10; //??
        let N = self.dim.1.row;
        match N {
            1 | 2 | 3 => self.inv(),
            _ => {
                let parts = Graph::<Cycle, Undirected>::from_matrix(self).tree_partition();
                if parts.len() < SPARSITY_THRESHHOLD { return self.invert_by_block() }
                //TODO may need case for block diagonal matrices

                m
            }
        }
    }
}

pub trait Ring: Mul<Output = Self> + Div<Output = Self> + Sub<Output = Self> + Add<Output = Self> + Zero + One + Neg<Output = Self> {}

pub trait Determinant<T> {
    fn det(&self) -> T;
}

pub trait Inverse<T> : Determinant<T> {
    fn inv(&self) -> Matrix<T>;
}

impl<T: Copy, const N: usize, const M: usize> From<[[T; M]; N]> for Matrix<T> {
    fn from(value: [[T; M]; N]) -> Self {
        Self{m: value.into_iter().flatten().collect::<Vec<T>>().into_boxed_slice(), dim: Dimension{row: N, col: M}}
    }
}

impl<T> From<Vec<Vec<T>>> for Matrix<T> {
    fn from(value: Vec<Vec<T>>) -> Self {
        let dim = Dimension{row: value.len(), col: value[0].len()};
        Matrix{m: value.into_iter().flatten().collect::<Vec<T>>().into_boxed_slice(), dim}
    }
}

impl<'a, 'b: 'a, T> From<&'b Matrix<T>> for SubMatrix<'a, T> {
    fn from(value: &'b Matrix<T>) -> Self {
        SubMatrix { m: &value, dim: ((0,0), value.dim) }
    }
}

impl<'a, T: Copy> From<SubMatrix<'a, T>> for Matrix<T> {
    fn from(value: SubMatrix<'a, T>) -> Self {
        unsafe {
            let mut tmp = Matrix::new_uninit(value.dim.1);
            tmp.copy_submatrix(value, (0,0));
            tmp
        }
    }
}

impl<'a, 'b, T: Mul<Output = T> + Add<Output = T> + Copy> Mul<SubMatrix<'a, T>> for SubMatrix<'b, T>{
    type Output = Matrix<T>;

    // TODO: sparse matrix multiplication?
    fn mul(self, rhs: SubMatrix<'a, T>) -> Self::Output {
        assert!(self.dim.1.col == rhs.dim.1.row);
        let mut m = Vec::with_capacity(self.dim.1.row * rhs.dim.1.col);
        for i in 0..self.dim.1.row {
            for j in 0..rhs.dim.1.col {
                let mut sum = self.get((i, 0)) * rhs.get((0, j));
                for k in 1..self.dim.1.col {
                    sum = sum + self.get((i, k)) * rhs.get((k, j));
                }
                m.push(sum)
            }
        }
        Matrix{m: m.into_boxed_slice(), dim: Dimension{row: self.dim.1.row, col: rhs.dim.1.col}}
    }

}

impl<T: Mul<Output = T> + Add<Output = T> + Copy> Mul<Matrix<T>> for Matrix<T>{
    type Output = Matrix<T>;

    fn mul(self, rhs: Matrix<T>) -> Self::Output {
        Into::<SubMatrix<'_, T>>::into(&self) * Into::<SubMatrix<'_, T>>::into(&rhs)
    }
}

impl<'a, T: Mul<Output = T> + Add<Output = T> + Copy> Mul<Matrix<T>> for SubMatrix<'a, T> {
    type Output = Matrix<T>;
    fn mul(self, rhs: Matrix<T>) -> Self::Output {
        self * Into::<SubMatrix<'_, T>>::into(&rhs)
    }
}

impl<'a, T: Mul<Output = T> + Add<Output = T> + Copy> Mul<SubMatrix<'a, T>> for Matrix<T> {
    type Output = Matrix<T>;
    fn mul(self, rhs: SubMatrix<'a, T>) -> Self::Output {
        Into::<SubMatrix<'_, T>>::into(&self) * rhs
    }
}

impl<'a, 'b, T: Add<Output = T> + Copy> Add<SubMatrix<'a, T>> for SubMatrix<'b, T> {
    type Output = Matrix<T>;

    fn add(self, rhs: SubMatrix<'a, T>) -> Self::Output {
        assert!(self.dim.1.col == rhs.dim.1.col && self.dim.1.row == rhs.dim.1.col);
        let mut m = Vec::with_capacity(self.dim.1.row * rhs.dim.1.col);
        for i in 0..self.dim.1.row {
            for j in 0..rhs.dim.1.col {
                m.push(self.get((i,j)) + rhs.get((i,j)));
            }
        }
        Matrix{m: m.into_boxed_slice(), dim: Dimension{row: self.dim.1.row, col: self.dim.1.col}}
    }
}

impl<'a, 'b, T: Sub<Output = T> + Copy> Sub<SubMatrix<'a, T>> for SubMatrix<'b, T> {
    type Output = Matrix<T>;

    fn sub(self, rhs: SubMatrix<'a, T>) -> Self::Output {
        assert!(self.dim.1.col == rhs.dim.1.col && self.dim.1.row == rhs.dim.1.col);
        let mut m = Vec::with_capacity(self.dim.1.row * rhs.dim.1.col);
        for i in 0..self.dim.1.row {
            for j in 0..rhs.dim.1.col {
                m.push(self.get((i,j)) - rhs.get((i,j)));
            }
        }
        Matrix{m: m.into_boxed_slice(), dim: Dimension{row: self.dim.1.row, col: self.dim.1.col}}
    }
}

impl<T: Add<Output = T> + Copy> Add<Matrix<T>> for Matrix<T>{
    type Output = Matrix<T>;

    fn add(self, rhs: Matrix<T>) -> Self::Output {
        Into::<SubMatrix<'_, T>>::into(&self) + Into::<SubMatrix<'_, T>>::into(&rhs)
    }
}

impl<T: Sub<Output = T> + Copy> Sub<Matrix<T>> for Matrix<T>{
    type Output = Matrix<T>;

    fn sub(self, rhs: Matrix<T>) -> Self::Output {
        Into::<SubMatrix<'_, T>>::into(&self) - Into::<SubMatrix<'_, T>>::into(&rhs)
    }
}

impl<'a, T: Add<Output = T> + Copy> Add<Matrix<T>> for SubMatrix<'a, T> {
    type Output = Matrix<T>;
    fn add(self, rhs: Matrix<T>) -> Self::Output {
        self + Into::<SubMatrix<'_, T>>::into(&rhs)
    }
}

impl<'a, T: Add<Output = T> + Copy> Add<SubMatrix<'a, T>> for Matrix<T> {
    type Output = Matrix<T>;
    fn add(self, rhs: SubMatrix<'a, T>) -> Self::Output {
        Into::<SubMatrix<'_, T>>::into(&self) + rhs
    }
}

impl<'a, T: Sub<Output = T> + Copy> Sub<SubMatrix<'a, T>> for Matrix<T> {
    type Output = Matrix<T>;
    fn sub(self, rhs: SubMatrix<'a, T>) -> Self::Output {
        Into::<SubMatrix<'_, T>>::into(&self) - rhs
    }
}

impl<'a, T: Sub<Output = T> + Copy> Sub<Matrix<T>> for SubMatrix<'a, T> {
    type Output = Matrix<T>;
    fn sub(self, rhs: Matrix<T>) -> Self::Output {
        self - Into::<SubMatrix<'_, T>>::into(&rhs)
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

impl<T: Div<Output = T> + Copy> Div<T> for Matrix<T> {
    type Output = Matrix<T>;
    fn div(self, rhs: T) -> Self::Output {
        Into::<SubMatrix<'_, T>>::into(&self) / rhs
    }
}

impl<T: Mul<Output = T> + Copy> Mul<T> for Matrix<T> {
    type Output = Matrix<T>;
    fn mul(self, rhs: T) -> Self::Output {
        Into::<SubMatrix<'_, T>>::into(&self) * rhs
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
        self.m.iter().zip(other.m.iter()).find(|(a, b)| **a != **b).is_none()
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
            2 => Into::<Matrix<T>>::into([[self.get((1,1)), -self.get((0,1))], [-self.get((1,0)), self.get((0,0))]]) / self.det(),
            3 => self.invert_3x3(),
            _ => self.invert_by_block(),
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
    use std::{io::{BufRead, Read}, vec};

    use super::*;
    impl Ring for i32 {}
    impl Ring for f64 {}

    const FLOAT_COMPARE: f64 = 0.00001;

    fn load_matrix<T: FromStr>(file: String) -> Matrix<T> {
        let mut f = std::fs::File::open("matrices/".to_owned() + &file + ".mat").expect(("failed to open matrix: ".to_owned() + &file).as_str());
        let mut buf = vec![];
        f.read_to_end(&mut buf).expect("failed to read file");
        buf.split(|c| *c==b';')
            .filter_map(|line| if line.len() > 0 { Some(line.split(|c| *c==b',')
                .filter_map(|c| if c.len() > 0 {
                    Some(FromStr::from_str(&String::from_utf8_lossy(c).into_owned()).ok().unwrap())
                } else {None}).collect()) }
            else {None}).collect::<Vec<Vec<T>>>().into()
    }

    fn load_determinant<T: FromStr>(test: String) -> T {
        let mut f = std::fs::File::open("matrices/".to_owned() + &test + ".mat.inf").expect(("failed to open matrix: ".to_owned() + &test).as_str());
        let mut buf = vec![];
        f.read_to_end(&mut buf).expect("failed to read .inf");
        let mut lines = buf.lines();
        let det = lines.next().unwrap().expect("failed to read line");
        FromStr::from_str(&det).ok().unwrap()
    }

    #[test]
    fn test_eq() {
        let a: Matrix<i32> = [[1, -1], [1, 1]].into();
        let b = [[1, -1], [1, 1]].into();
        assert!(a == b);
    }

    #[test]
    fn test_mul() {
        let a: Matrix<i32> = [[1_i32, -1], [1, 1]].into();
        let b: Matrix<i32> = [[-1, -1], [2, 1]].into();
        let c = [[-3, -2], [1, 0]].into();
        let d: Matrix<i32> = [[3,4,5], [-1,0,1]].into();
        let e = [[4, 4, 4], [2, 4, 6]].into();
        assert!(a.clone() * b == c);
        assert!(a * d == e);
    }

    #[test]
    fn test_add_sub() {
        let a: Matrix<i32> = [[1_i32, -1], [1, 1]].into();
        let b: Matrix<i32> = [[-1, -1], [2, 1]].into();
        let c = [[0, -2], [3, 2]].into();
        let d = [[2, 0], [-1, 0]].into();
        assert_eq!(a.clone() + b.clone(), c);
        assert_eq!(a.clone() - b.clone(), d);
        assert_eq!(a + (b * -1), d);
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
    fn test_det_large() {
        let a: Matrix<i32> = [[0, 0, -3, 0, 0, 0, -3, 0, 4, 0], [0, 0, 0, 0, 0, -5, 6, -10, -2, 0], [-7, 0, 0, 0, 10, 0, 0, 1, 0, 0], [0, 8, 0, 0, 5, 0, -7, 0, 0, -3], [0, 0, 0, 0, 0, -3, -10, 4, 6, 0], [5, 0, -1, 0, 8, 1, 0, 0, -7, 0], [0, 3, 0, 0, 0, 0, 0, 2, 0, -1], [0, 0, 0, 5, -5, 0, 0, -9, 3, 0], [0, 0, 6, 0, 0, 0, 8, 0, 0, 5], [1, 0, 0, -5, 0, -7, 0, 1, 0, 0]].into();
        assert!(a.det() == 7396440)
    }

    #[test]
    fn test_det_float() {
        for name in ["test1", "test2", "test3", "test9", "test10"] {
            let a: Matrix<f64> = load_matrix(String::from(name));
            println!("{:?}", a.det());
            let expected = load_determinant::<f64>(String::from(name));
            assert!(expected - a.det() < FLOAT_COMPARE, "failed: {:?}", name);
        }
    }

    #[test]
    fn test_inv_base() {
        let a: Matrix<i32> = [[2, 1], [1, 1]].into();
        assert!(a.clone() * a.clone().inv() == [[1,0],[0,1]].into());
        assert!(a.clone().inv() * a == [[1,0],[0,1]].into());
        let a: Matrix<i32> = [[-2, -7], [1, 4]].into();
        assert!(a.clone() * a.clone().inv() == [[1,0],[0,1]].into());
        assert!(a.clone().inv() * a == [[1,0],[0,1]].into());
    }

    #[test]
    fn test_from_submatrices() {
        let a: Matrix<i32> = [[1,2],[3,4]].into();
        let b: Matrix<i32> = [[5,6], [7,8]].into();
        let c: Matrix<i32> = [[9,10],[11,12]].into();
        let d: Matrix<i32> = [[13,14],[15,16]].into();
        let a = Matrix::from_submatrices((&a).into(), (&b).into(), (&c).into(), (&d).into());
        let e: Matrix<i32> = SubMatrix::from_matrix(&a, ((1,1),Dimension{row: 2, col: 2})).into();
        assert_eq!(a, [[1, 2, 5, 6], [3, 4, 7, 8], [9, 10, 13, 14], [11, 12, 15, 16]].into());
        assert_eq!(e, [[4,7],[10,13]].into());
    }

    #[test]
    fn test_inv_med() {
        for name in ["test2", "test3", "test4", "test5", "test6", "test9"] {
            let a: Matrix<f64> = load_matrix(String::from(name));
            let ai = a.inv();
            println!("testing {:?}", name);
            assert!(SubMatrix::nearly_equal(&(&(a.clone() * ai.clone())).into(), &(&Matrix::identity(a.dim.col)).into()));
            assert!(SubMatrix::nearly_equal(&(&(ai * a.clone())).into(), &(&Matrix::identity(a.dim.col)).into()));
        }
    }

    // #[test]
    // fn test_gray_code() {
    //     assert_eq!((0..4).map(|n| gray_code(n)).collect::<Vec<u32>>(), &[0,1,3,2])
    // }
}
