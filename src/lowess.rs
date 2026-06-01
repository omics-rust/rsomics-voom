//! Cleveland's LOWESS (locally weighted scatterplot smoothing), matching the
//! `lowess()` in R base `stats` — a port of Cleveland's `clowess`.
//! Reference: W. S. Cleveland (1979), JASA 74:829-836.
//!
//! `x` must be sorted ascending.

pub fn lowess(x: &[f64], y: &[f64], f: f64, iter: usize, delta: f64) -> Vec<f64> {
    let n = x.len();
    let mut ys = vec![0.0; n];
    if n < 2 {
        if n == 1 {
            ys[0] = y[0];
        }
        return ys;
    }

    let ns = ((f * n as f64 + 1e-7) as usize).clamp(2, n);

    let mut rw = vec![0.0; n];
    let mut res = vec![0.0; n];
    let mut w = vec![0.0; n];

    for it in 0..=iter {
        let mut nleft = 0usize;
        let mut nright = ns - 1;
        let mut last: isize = -1;
        let mut i = 0usize;

        loop {
            while nright < n - 1 {
                let d1 = x[i] - x[nleft];
                let d2 = x[nright + 1] - x[i];
                if d1 <= d2 {
                    break;
                }
                nleft += 1;
                nright += 1;
            }

            let robust = if it > 0 { Some(&rw[..]) } else { None };
            ys[i] = lowest(x, y, x[i], (nleft, nright), &mut w, robust);

            if last < i as isize - 1 {
                let denom = x[i] - x[last as usize];
                for j in (last + 1) as usize..i {
                    let alpha = (x[j] - x[last as usize]) / denom;
                    ys[j] = alpha * ys[i] + (1.0 - alpha) * ys[last as usize];
                }
            }

            last = i as isize;
            let cut = x[last as usize] + delta;
            i = (last + 1) as usize;
            while i < n {
                if x[i] > cut {
                    break;
                }
                if x[i] == x[last as usize] {
                    ys[i] = ys[last as usize];
                    last = i as isize;
                }
                i += 1;
            }
            i = ((last + 1) as usize).max(i.saturating_sub(1));

            if last >= (n - 1) as isize {
                break;
            }
        }

        if it == iter {
            break;
        }

        for k in 0..n {
            res[k] = y[k] - ys[k];
        }
        let mut abs_res: Vec<f64> = res.iter().map(|r| r.abs()).collect();
        abs_res.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let m1 = n / 2;
        let cmad = if n.is_multiple_of(2) {
            3.0 * (abs_res[m1 - 1] + abs_res[m1])
        } else {
            6.0 * abs_res[m1]
        };

        let sum_abs: f64 = abs_res.iter().sum();
        if cmad < 1e-7 * sum_abs {
            break;
        }

        let c9 = 0.999 * cmad;
        let c1 = 0.001 * cmad;
        for k in 0..n {
            let r = res[k].abs();
            if r <= c1 {
                rw[k] = 1.0;
            } else if r <= c9 {
                let t = r / cmad;
                let u = 1.0 - t * t;
                rw[k] = u * u;
            } else {
                rw[k] = 0.0;
            }
        }
    }

    ys
}

fn lowest(
    x: &[f64],
    y: &[f64],
    xs: f64,
    window: (usize, usize),
    w: &mut [f64],
    rw: Option<&[f64]>,
) -> f64 {
    let (nleft, nright) = window;
    let n = x.len();
    let range = x[n - 1] - x[0];
    let h = (xs - x[nleft]).max(x[nright] - xs);
    let h9 = 0.999 * h;
    let h1 = 0.001 * h;

    let mut a = 0.0;
    for j in nleft..=nright {
        w[j] = 0.0;
        let r = (x[j] - xs).abs();
        if r <= h9 {
            if r > h1 {
                let t = r / h;
                let cube = 1.0 - t * t * t;
                w[j] = cube * cube * cube; // tricube
            } else {
                w[j] = 1.0;
            }
            if let Some(rw) = rw {
                w[j] *= rw[j];
            }
            a += w[j];
        }
    }

    if a <= 0.0 {
        return y[nleft];
    }

    for wk in w.iter_mut().take(nright + 1).skip(nleft) {
        *wk /= a;
    }

    if h > 0.0 {
        let mut a_mean = 0.0;
        for k in nleft..=nright {
            a_mean += w[k] * x[k];
        }
        let b = xs - a_mean;
        let mut c = 0.0;
        for k in nleft..=nright {
            c += w[k] * (x[k] - a_mean) * (x[k] - a_mean);
        }
        if c.sqrt() > 0.001 * range {
            let bb = b / c;
            for k in nleft..=nright {
                w[k] *= bb * (x[k] - a_mean) + 1.0;
            }
        }
    }

    let mut fit = 0.0;
    for k in nleft..=nright {
        fit += w[k] * y[k];
    }
    fit
}

/// `approxfun(x, y, rule = 2, ties = list("ordered", mean))`: equal-x groups
/// collapsed to their mean y once, then strictly-increasing x for O(log n)
/// linear interpolation clamped to the endpoints.
pub struct Trend {
    x: Vec<f64>,
    y: Vec<f64>,
}

impl Trend {
    pub fn new(lx: &[f64], ly: &[f64]) -> Self {
        let mut x = Vec::new();
        let mut y = Vec::new();
        let mut i = 0;
        while i < lx.len() {
            let xv = lx[i];
            let mut j = i;
            let mut sum = 0.0;
            while j < lx.len() && lx[j] == xv {
                sum += ly[j];
                j += 1;
            }
            x.push(xv);
            y.push(sum / (j - i) as f64);
            i = j;
        }
        Trend { x, y }
    }

    pub fn eval(&self, xout: f64) -> f64 {
        let last = self.x.len() - 1;
        if xout <= self.x[0] {
            return self.y[0];
        }
        if xout >= self.x[last] {
            return self.y[last];
        }
        let hi = self.x.partition_point(|&v| v <= xout);
        let lo = hi - 1;
        let (x0, x1) = (self.x[lo], self.x[hi]);
        self.y[lo] + (self.y[hi] - self.y[lo]) * (xout - x0) / (x1 - x0)
    }
}
