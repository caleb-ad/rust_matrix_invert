import random
from numpy import matrix
from numpy.linalg import det


def gen_inverse_tests():
    for (test, size, sparsity) in [("test1", 3, .1), ("test2", 5, .01), ("test3", 9, .75), ("test4", 10, .5)]:
        m1 = gen_matrix(size, sparsity=sparsity)
        write_rust_matrix("matrices/" + test + ".mat", m1)
        write_info("matrices/" + test + ".mat.inf", m1)
        # m1_inv = matrix(m1).I
        # write_rust_matrix("matrices/" + test + ".mat.inv", m1_inv.tolist())


def gen_matrix(n, sparsity=.5):
    return [[random.randint(-10, 10) if random.randint(0,100) > 100*sparsity else 0 for x in range(n)] for y in range(n)]

def write_rust_matrix(file_name, m):
    with open(file_name, mode='w') as f:
        for row in m:
            for entry in row:
                f.write("{:.10f},".format(entry))
            f.write(";")

def write_info(file_name, m):
    with open(file_name, mode='w') as f:
        f.write("{:.10f}\n".format(det(m)))


if __name__ == "__main__":
    gen_inverse_tests()
