def remap_colors(old_delta):
    if old_delta <= 180:
        new_delta = old_delta
    else:
        new_delta = 360 - old_delta
    if new_delta == 0:
        new_delta = 0.0001
    else:
        new_delta = new_delta/2.0
    return new_delta
