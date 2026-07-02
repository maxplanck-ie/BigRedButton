def edit_distance(a: str, b: str, max_dist: int | None = None) -> int:
    if a == b:
        return 0

    la = len(a)
    lb = len(b)
    if max_dist is not None and abs(la - lb) > max_dist:
        return max_dist + 1

    if la < lb:
        a, b = b, a
        la, lb = lb, la

    previous = list(range(lb + 1))
    for i in range(1, la + 1):
        current = [i] + [0] * lb
        best = current[0]
        ai = a[i - 1]
        for j in range(1, lb + 1):
            cost = 0 if ai == b[j - 1] else 1
            current[j] = min(
                previous[j] + 1,
                current[j - 1] + 1,
                previous[j - 1] + cost,
            )
            if current[j] < best:
                best = current[j]
        if max_dist is not None and best > max_dist:
            return max_dist + 1
        previous = current

    return previous[lb]


def eval(a: str, b: str, max_dist: int | None = None) -> int:
    return edit_distance(a, b, max_dist=max_dist)
