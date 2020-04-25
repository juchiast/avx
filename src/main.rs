#[cfg(target_arch = "x86")]
use std::arch::x86::*;
#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::*;

fn factor(n: u32) -> (u32, u32) {
    let s = n.trailing_zeros();
    (s, n >> s)
}

unsafe fn x2mod(a: __m256i, modulo: __m256i) -> __m256i {
    // x = a + a
    // return x - (x > N? 1 : 0) * N
    let x = _mm256_add_epi64(a, a);
    _mm256_sub_epi64(x, _mm256_and_si256(_mm256_cmpgt_epi64(x, modulo), modulo))
}

unsafe fn addmod(a: __m256i, b: __m256i, modulo: __m256i) -> __m256i {
    let x = _mm256_add_epi64(a, b);
    _mm256_sub_epi64(x, _mm256_and_si256(_mm256_cmpgt_epi64(x, modulo), modulo))
}

unsafe fn mulmod(mut a: __m256i, mut b: __m256i, modulo: __m256i) -> __m256i {
    let mut result = _mm256_set1_epi64x(0);
    let one = _mm256_set1_epi64x(1);
    for _ in 0..32 {
        result = addmod(
            result,
            _mm256_mul_epu32(a, _mm256_and_si256(b, one)),
            modulo,
        );
        a = x2mod(a, modulo);
        b = _mm256_srli_epi64(b, 1)
    }
    result
}

unsafe fn pow(mut a: __m256i, mut d: u32, modulo: __m256i) -> __m256i {
    let mut result = _mm256_set1_epi64x(1);
    while d > 0 {
        if d & 1 == 1 {
            result = mulmod(result, a, modulo);
        }
        a = mulmod(a, a, modulo);
        d >>= 1;
    }
    result
}

unsafe fn witness_test(n: i64, s: u32, d: u32) -> bool {
    let witness = _mm256_set_epi64x(2 % n, 7 % n, 61 % n, 0);
    let modulo = _mm256_set1_epi64x(n);
    let mut p = pow(witness, d, modulo);
    let mut b = _mm256_cmpeq_epi64(p, _mm256_set_epi64x(1, 1, 1, 0));
    let temp = _mm256_set1_epi64x(n - 1);
    for _ in 0..s {
        b = _mm256_or_si256(b, _mm256_cmpeq_epi64(p, temp));
        p = mulmod(p, p, modulo);
    }
    _mm256_movemask_epi8(b) == -1
}

fn miller(n: u32) -> bool {
    // if n < 2 { return false; }
    // if n % 2 == 0 { return false; }
    if n == 2 || n == 7 || n == 61 { return true; }
    let (s, d) = factor(n - 1);
    unsafe { witness_test(n as i64, s, d) }
}

fn main() {
    let to = (1i64 << 32) - 1;
    for n in (3i64..to).step_by(2) {
        if miller_rabin::is_prime(&n, 100) != miller(n as u32) {
            println!("{}", n);
        }
    }
    /*
    use std::io::Read;
    let mut input = String::new();
    std::io::stdin().read_to_string(&mut input).ok();
    let mut tests: Vec<(u64, usize)> = input
        .split_whitespace()
        .skip(1)
        .enumerate()
        .map(|(i, x)| (x.parse().unwrap(), i))
        .collect();
    tests.sort_unstable();
    let mut ans: Vec<u32> = Vec::new();
    ans.resize(tests.len(), 0);

    let mut last_prime = 2u64;
    let mut bound = last_prime + 1;
    for (n, i) in &tests {
        let mut p = if n&1 == 0 { n - 1 } else { n - 2 };
        while p >= bound {
            if miller(p as u32) {
                last_prime = p;
                break;
            }
            p -= 2;
        }

        ans[*i] = last_prime as u32;
        bound = *n;
    }


    let mut output = String::new();
    for x in ans {
        output.push_str(&x.to_string());
        output.push('\n');
    }
    print!("{}", output);
    */
}
