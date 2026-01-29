from common import *

@dataclass(frozen=True)
class AffineFunc:
    a: int
    b: int
    def __post_init__(self):
        if math.gcd(self.a, P) != 1:
            raise ValueError(f"a ({self.a}) must be coprime with P ({P})")

    def f(self, x: int) -> int:
        return (self.a * x + self.b) % P

    def inv(self):
        try:
            a_inv = pow(self.a, -1, P)
        except ValueError:
             raise ValueError("Inverse does not exist")
        b_new = (-1 * a_inv * self.b) % P
        return AffineFunc(a_inv, b_new)

def composite_affine(left, right):
    a_new = (left.a * right.a) % P
    b_new = (left.a * right.b + left.b) % P
    return AffineFunc(a_new, b_new)

def generate_random_affine_function():
    prime_list = [i for i in range(P) if math.gcd(i,P) == 1]
    a = random.choice(prime_list)
    b = random.choice(range(P))
    return AffineFunc(a, b)

def commute_with(f):
    a_1 = f.a
    b_1 = f.b
    candidates = []
    for a_2 in range(0,P):
        if math.gcd(a_2,P)!=1:
            continue
        A = a_1-1
        C = ((a_2-1)*b_1)%P
        d = math.gcd(A, P)
        if a_2==1:
            if a_1==1:
                for i in range(P):
                    b_2 = i
                    candidates.append(AffineFunc(a_2, b_2))
                continue
            else:
                if d==1:
                    b_2 = 0
                    candidates.append(AffineFunc(a_2, b_2))
                    continue
                else:
                    P_prime = int(P/d)
                    for i in range(d):
                        b_2 = i*P_prime
                        candidates.append(AffineFunc(a_2, b_2))
                    continue
        else:
            if a_1==1:
                if C==0:
                    for i in range(P):
                        b_2 = i
                        candidates.append(AffineFunc(a_2, b_2))
                    continue
            else:
                if d==1:
                    b_2 = (pow(A, -1, P) * C)%P
                    candidates.append(AffineFunc(a_2, b_2))
                    continue
                else:
                    if C%d==0:
                        A_prime = int(A/d)
                        P_prime = int(P/d)
                        C_prime = int(C/d)
                        x_0 = (C_prime*pow(A_prime, -1, P_prime))%P_prime
                        for i in range(d):
                            b_2 = x_0+i*P_prime
                            candidates.append(AffineFunc(a_2, b_2))
                        continue
    return candidates