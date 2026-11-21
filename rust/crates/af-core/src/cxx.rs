//! Behaviour of the C++ standard library that results depend on. `sort_by`
//! is libstdc++'s `std::sort` (introsort: median-of-three quicksort to depth
//! 2 log2 n, heap sort beyond, then insertion sort over runs of 16), so equal
//! elements end in the same order as in the C++ build, where an algorithm
//! breaks ties by position (neighbour frames, region seeds).

const THRESHOLD: usize = 16;

/// `std::sort(v.begin(), v.end(), less)` with libstdc++'s exact element order.
pub fn sort_by<T, F: FnMut(&T, &T) -> bool>(v: &mut [T], mut less: F) {
    let n = v.len();
    if n == 0 {
        return;
    }
    let depth = 2 * (usize::BITS - 1 - n.leading_zeros()) as usize;
    introsort_loop(v, 0, n, depth, &mut less);
    final_insertion_sort(v, 0, n, &mut less);
}

fn introsort_loop<T, F: FnMut(&T, &T) -> bool>(v: &mut [T], first: usize, mut last: usize, mut depth: usize, less: &mut F) {
    while last - first > THRESHOLD {
        if depth == 0 {
            partial_sort_whole(v, first, last, less);
            return;
        }
        depth -= 1;
        let cut = unguarded_partition_pivot(v, first, last, less);
        introsort_loop(v, cut, last, depth, less);
        last = cut;
    }
}

fn unguarded_partition_pivot<T, F: FnMut(&T, &T) -> bool>(v: &mut [T], first: usize, last: usize, less: &mut F) -> usize {
    let mid = first + (last - first) / 2;
    move_median_to_first(v, first, first + 1, mid, last - 1, less);
    unguarded_partition(v, first + 1, last, first, less)
}

fn move_median_to_first<T, F: FnMut(&T, &T) -> bool>(v: &mut [T], result: usize, a: usize, b: usize, c: usize, less: &mut F) {
    if less(&v[a], &v[b]) {
        if less(&v[b], &v[c]) {
            v.swap(result, b);
        } else if less(&v[a], &v[c]) {
            v.swap(result, c);
        } else {
            v.swap(result, a);
        }
    } else if less(&v[a], &v[c]) {
        v.swap(result, a);
    } else if less(&v[b], &v[c]) {
        v.swap(result, c);
    } else {
        v.swap(result, b);
    }
}

fn unguarded_partition<T, F: FnMut(&T, &T) -> bool>(v: &mut [T], mut first: usize, mut last: usize, pivot: usize, less: &mut F) -> usize {
    loop {
        while less(&v[first], &v[pivot]) {
            first += 1;
        }
        last -= 1;
        while less(&v[pivot], &v[last]) {
            last -= 1;
        }
        if first >= last {
            return first;
        }
        v.swap(first, last);
        first += 1;
    }
}

fn final_insertion_sort<T, F: FnMut(&T, &T) -> bool>(v: &mut [T], first: usize, last: usize, less: &mut F) {
    if last - first > THRESHOLD {
        insertion_sort(v, first, first + THRESHOLD, less);
        for i in first + THRESHOLD..last {
            unguarded_linear_insert(v, i, less);
        }
    } else {
        insertion_sort(v, first, last, less);
    }
}

fn insertion_sort<T, F: FnMut(&T, &T) -> bool>(v: &mut [T], first: usize, last: usize, less: &mut F) {
    if first == last {
        return;
    }
    for i in first + 1..last {
        if less(&v[i], &v[first]) {
            // Move the element to the front, shifting the rest up by one.
            v[first..=i].rotate_right(1);
        } else {
            unguarded_linear_insert(v, i, less);
        }
    }
}

fn unguarded_linear_insert<T, F: FnMut(&T, &T) -> bool>(v: &mut [T], mut last: usize, less: &mut F) {
    // Shift the element at `last` down while it is less than its predecessor.
    while less(&v[last], &v[last - 1]) {
        v.swap(last, last - 1);
        last -= 1;
    }
}

// Heap sort of [first, last) as std::__partial_sort(first, last, last):
// __heap_select (make_heap) followed by __sort_heap.
fn partial_sort_whole<T, F: FnMut(&T, &T) -> bool>(v: &mut [T], first: usize, last: usize, less: &mut F) {
    let s = &mut v[first..last];
    make_heap(s, less);
    // __heap_select over [middle, last) is empty when middle == last.
    let mut len = s.len();
    while len > 1 {
        len -= 1;
        // __pop_heap(first, last, last): value = *last; *last = *first; adjust_heap(first, 0, len, value).
        s.swap(0, len);
        adjust_heap_with_hole(&mut s[..len], 0, less);
    }
}

fn make_heap<T, F: FnMut(&T, &T) -> bool>(s: &mut [T], less: &mut F) {
    let len = s.len();
    if len < 2 {
        return;
    }
    let mut parent = (len - 2) / 2;
    loop {
        adjust_heap_with_hole(s, parent, less);
        if parent == 0 {
            return;
        }
        parent -= 1;
    }
}

// std::__adjust_heap(first, hole, len, value) with value = s[hole]: sift the
// hole down along the larger children, then __push_heap the value back up.
fn adjust_heap_with_hole<T, F: FnMut(&T, &T) -> bool>(s: &mut [T], hole: usize, less: &mut F) {
    let len = s.len();
    let top = hole;
    let mut hole = hole;
    let mut child = hole;
    // The value travels with swaps, which reproduces the hole moves exactly.
    while len > 0 && child < (len - 1) / 2 {
        child = 2 * (child + 1);
        if less(&s[child], &s[child - 1]) {
            child -= 1;
        }
        s.swap(hole, child);
        hole = child;
    }
    if len % 2 == 0 && len >= 2 && child == (len - 2) / 2 {
        child = 2 * (child + 1);
        s.swap(hole, child - 1);
        hole = child - 1;
    }
    // __push_heap(first, hole, top, value)
    while hole > top {
        let parent = (hole - 1) / 2;
        if !less(&s[parent], &s[hole]) {
            break;
        }
        s.swap(hole, parent);
        hole = parent;
    }
}

/// `std::nth_element(v.begin(), v.begin() + nth, v.end(), less)` with
/// libstdc++'s exact element order (introselect: median-of-three
/// partitioning to depth 2 log2 n, then heap selection; insertion sort for
/// the last three or fewer).
pub fn nth_element_by<T, F: FnMut(&T, &T) -> bool>(v: &mut [T], nth: usize, mut less: F) {
    let n = v.len();
    if n == 0 || nth >= n {
        return;
    }
    let mut depth = 2 * (usize::BITS - 1 - n.leading_zeros()) as usize;
    let (mut first, mut last) = (0, n);
    while last - first > 3 {
        if depth == 0 {
            heap_select(v, first, nth + 1, last, &mut less);
            // Place the nth largest element in its final position.
            v.swap(first, nth);
            return;
        }
        depth -= 1;
        let cut = unguarded_partition_pivot(v, first, last, &mut less);
        if cut <= nth {
            first = cut;
        } else {
            last = cut;
        }
    }
    insertion_sort(v, first, last, &mut less);
}

/// `std::partial_sort(v.begin(), v.begin() + middle, v.end(), less)` with
/// libstdc++'s exact element order: heap selection, then heap sort of the
/// first `middle` elements.
pub fn partial_sort_by<T, F: FnMut(&T, &T) -> bool>(v: &mut [T], middle: usize, mut less: F) {
    let middle = middle.min(v.len());
    // With middle == first libstdc++ still runs the selection loop (it swaps
    // smaller elements to the front).
    heap_select(v, 0, middle, v.len(), &mut less);
    let s = &mut v[..middle];
    let mut len = s.len();
    while len > 1 {
        len -= 1;
        s.swap(0, len);
        adjust_heap_with_hole(&mut s[..len], 0, &mut less);
    }
}

// std::__heap_select(first, middle, last): a heap of [first, middle), then
// every smaller element of [middle, last) replaces its top.
fn heap_select<T, F: FnMut(&T, &T) -> bool>(v: &mut [T], first: usize, middle: usize, last: usize, less: &mut F) {
    make_heap(&mut v[first..middle], less);
    for i in middle..last {
        if less(&v[i], &v[first]) {
            // __pop_heap(first, middle, i): the value of *i sifts down from the top.
            v.swap(first, i);
            adjust_heap_with_hole(&mut v[first..middle], 0, less);
        }
    }
}

/// std::stable_sort as libstdc++ (GCC 15) does it with its temporary buffer:
/// identical results even for comparators that are not strict weak orders
/// (tolerance comparisons).
pub fn stable_sort_by<T: Clone, F: FnMut(&T, &T) -> bool>(v: &mut [T], mut less: F) {
    let len = v.len();
    if len == 0 {
        return;
    }
    let half = len.div_ceil(2);
    let mut buffer: Vec<T> = v[..half].to_vec();
    // __stable_sort_adaptive(first, first + half, last, buffer)
    merge_sort_with_buffer(v, 0, half, &mut buffer, &mut less);
    merge_sort_with_buffer(v, half, len, &mut buffer, &mut less);
    merge_adaptive(v, 0, half, len, &mut buffer, &mut less);
}

const CHUNK: usize = 7;

fn merge_sort_with_buffer<T: Clone, F: FnMut(&T, &T) -> bool>(v: &mut [T], first: usize, last: usize, buffer: &mut Vec<T>, less: &mut F) {
    let len = last - first;
    let mut step = CHUNK;
    // __chunk_insertion_sort
    let mut f = first;
    while last - f >= step {
        insertion_sort(v, f, f + step, less);
        f += step;
    }
    insertion_sort(v, f, last, less);
    while step < len {
        {
            let src = v[first..last].to_vec();
            merge_sort_loop(&src, &mut buffer[..len], step, less);
        }
        step *= 2;
        {
            let src = buffer[..len].to_vec();
            merge_sort_loop(&src, &mut v[first..last], step, less);
        }
        step *= 2;
    }
}

/// __merge_sort_loop: runs of `step` from src merged pairwise into dst.
fn merge_sort_loop<T: Clone, F: FnMut(&T, &T) -> bool>(src: &[T], dst: &mut [T], step: usize, less: &mut F) {
    let two = 2 * step;
    let (mut first, last) = (0usize, src.len());
    let mut out = 0usize;
    while last - first >= two {
        out = move_merge(src, first, first + step, first + step, first + two, dst, out, less);
        first += two;
    }
    let step = (last - first).min(step);
    move_merge(src, first, first + step, first + step, last, dst, out, less);
}

#[allow(clippy::too_many_arguments)]
fn move_merge<T: Clone, F: FnMut(&T, &T) -> bool>(
    src: &[T],
    mut f1: usize,
    l1: usize,
    mut f2: usize,
    l2: usize,
    dst: &mut [T],
    mut out: usize,
    less: &mut F,
) -> usize {
    while f1 != l1 && f2 != l2 {
        if less(&src[f2], &src[f1]) {
            dst[out] = src[f2].clone();
            f2 += 1;
        } else {
            dst[out] = src[f1].clone();
            f1 += 1;
        }
        out += 1;
    }
    for item in &src[f1..l1] {
        dst[out] = item.clone();
        out += 1;
    }
    for item in &src[f2..l2] {
        dst[out] = item.clone();
        out += 1;
    }
    out
}

fn merge_adaptive<T: Clone, F: FnMut(&T, &T) -> bool>(v: &mut [T], first: usize, middle: usize, last: usize, buffer: &mut [T], less: &mut F) {
    let (len1, len2) = (middle - first, last - middle);
    if len1 <= len2 {
        // __move_merge_adaptive(buffer[..len1], v[middle..last]) into v[first..]
        buffer[..len1].clone_from_slice(&v[first..middle]);
        let (mut b, mut f2, mut out) = (0usize, middle, first);
        while b != len1 && f2 != last {
            if less(&v[f2], &buffer[b]) {
                v[out] = v[f2].clone();
                f2 += 1;
            } else {
                v[out] = buffer[b].clone();
                b += 1;
            }
            out += 1;
        }
        while b != len1 {
            v[out] = buffer[b].clone();
            b += 1;
            out += 1;
        }
    } else {
        // __move_merge_adaptive_backward(v[first..middle], buffer[..len2]) ending at last
        buffer[..len2].clone_from_slice(&v[middle..last]);
        if first == middle {
            v[first..last].clone_from_slice(&buffer[..len2]);
            return;
        }
        if len2 == 0 {
            return;
        }
        let (mut l1, mut l2, mut result) = (middle - 1, len2 - 1, last);
        loop {
            if less(&buffer[l2], &v[l1]) {
                result -= 1;
                v[result] = v[l1].clone();
                if l1 == first {
                    // move_backward(buffer[0..=l2]) ending at result
                    let n = l2 + 1;
                    v[result - n..result].clone_from_slice(&buffer[..n]);
                    return;
                }
                l1 -= 1;
            } else {
                result -= 1;
                v[result] = buffer[l2].clone();
                if l2 == 0 {
                    return;
                }
                l2 -= 1;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::{nth_element_by, sort_by, stable_sort_by};

    #[test]
    fn stable_sorts() {
        for n in [0usize, 1, 2, 6, 7, 8, 13, 14, 15, 28, 29, 57, 100, 1000] {
            let v: Vec<(i32, usize)> = (0..n).map(|i| (((i * 7919) % 13) as i32, i)).collect();
            let mut a = v.clone();
            stable_sort_by(&mut a, |x, y| x.0 < y.0);
            let mut b = v.clone();
            b.sort_by_key(|x| x.0);
            assert_eq!(a, b, "n = {n}");
        }
    }

    #[test]
    fn sorts() {
        let mut v: Vec<i32> = (0..1000).map(|i| (i * 7919) % 211).collect();
        sort_by(&mut v, |a, b| a < b);
        assert!(v.windows(2).all(|w| w[0] <= w[1]));
    }

    #[test]
    fn selects() {
        for n in [1usize, 2, 3, 4, 5, 17, 100, 1000] {
            for nth in [0, n / 3, n / 2, n - 1] {
                let mut v: Vec<i32> = (0..n as i32).map(|i| (i * 7919) % 97).collect();
                nth_element_by(&mut v, nth, |a, b| a < b);
                let pivot = v[nth];
                assert!(v[..nth].iter().all(|&x| x <= pivot) && v[nth..].iter().all(|&x| x >= pivot));
                let mut sorted = v.clone();
                sorted.sort();
                assert_eq!(pivot, sorted[nth]);
            }
        }
    }
}
