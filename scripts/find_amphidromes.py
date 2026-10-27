"""
Find every amphidromic point in NASA's GOT4.10c global ocean tide model.

Method: the tide at each grid point is A*cos(w t - G), with A = amplitude and
G = Greenwich phase lag. Walk counter-clockwise around every 0.5 deg grid cell
and add up the change in G. If G winds by +360 or -360 deg, the cell holds an
amphidromic point (+1 = high water circles counter-clockwise seen from above,
-1 = clockwise). The point is then located inside the cell by finding the zero
of the bilinear-interpolated complex tide Z = A*exp(-iG).

Usage:
    python -I find_amphidromes.py <GOT4.10c/grids_oceantide dir> <output .js>

Data: R. Ray, GOT4.10c, NASA GSFC (earth.gsfc.nasa.gov/geo/data/ocean-tide-models)
"""
import sys, json, base64
import numpy as np

CONSTITUENTS = ["m2", "s2", "k1", "o1"]
PERIOD_H = {"m2": 12.4206012, "s2": 12.0, "k1": 23.9344696, "o1": 25.8193417}
R_EARTH = 6371.0


def read_got_ascii(path):
    """Two grids (amplitude cm, phase deg) in GOT ascii, rows = latitude, 11F7.2."""
    with open(path) as f:
        lines = f.read().splitlines()
    grids, k = [], 0
    for _ in range(2):
        nlat, nlon = map(int, lines[k + 2].split())
        lat0, lat1 = map(float, lines[k + 3].split())
        lon0, lon1 = map(float, lines[k + 4].split())
        fill = float(lines[k + 5].split()[0])
        k += 7
        rows_per_lat = -(-nlon // 11)
        g = np.empty((nlat, nlon))
        for j in range(nlat):
            vals = []
            for ln in lines[k:k + rows_per_lat]:
                vals += [float(ln[c:c + 7]) for c in range(0, len(ln), 7) if ln[c:c + 7].strip()]
            g[j] = vals[:nlon]
            k += rows_per_lat
        grids.append(g)
    lat = np.linspace(lat0, lat1, nlat)
    lon = np.linspace(lon0, lon1, nlon)
    return lat, lon, grids[0], grids[1], fill


def wrap(d):
    return (d + 180.0) % 360.0 - 180.0


def find_points(lat, lon, amp, pha, fill):
    ocean = amp < fill - 1
    Z = np.where(ocean, amp * np.exp(-1j * np.radians(pha)), 0)
    nlat, nlon = amp.shape
    ip = (np.arange(nlon) + 1) % nlon
    # corners CCW in (east, north): (j,i) -> (j,i+1) -> (j+1,i+1) -> (j+1,i)
    c0, c1 = pha[:-1, :], pha[:-1, ip]
    c2, c3 = pha[1:, ip], pha[1:, :]
    ok = ocean[:-1, :] & ocean[:-1, ip] & ocean[1:, ip] & ocean[1:, :]
    wind = (wrap(c1 - c0) + wrap(c2 - c1) + wrap(c3 - c2) + wrap(c0 - c3)) / 360.0
    wind = np.where(ok, np.rint(wind), 0).astype(int)
    pts = []
    for j, i in zip(*np.nonzero(wind)):
        z00, z10, z11, z01 = Z[j, i], Z[j, ip[i]], Z[j + 1, ip[i]], Z[j + 1, i]
        u = v = 0.5  # Newton on bilinear Z(u,v) = 0
        for _ in range(30):
            zz = (1-u)*(1-v)*z00 + u*(1-v)*z10 + u*v*z11 + (1-u)*v*z01
            du = (1-v)*(z10 - z00) + v*(z11 - z01)
            dv = (1-u)*(z01 - z00) + u*(z11 - z10)
            J = np.array([[du.real, dv.real], [du.imag, dv.imag]])
            try:
                step = np.linalg.solve(J, [-zz.real, -zz.imag])
            except np.linalg.LinAlgError:
                break
            u, v = np.clip(u + step[0], 0, 1), np.clip(v + step[1], 0, 1)
        dlat, dlon = lat[1] - lat[0], lon[1] - lon[0]
        la = lat[j] + v * dlat
        lo = (lon[i] + u * dlon + 180) % 360 - 180
        pts.append({"lat": round(float(la), 3), "lon": round(float(lo), 3), "s": int(wind[j, i])})
    return pts, ocean


def add_land_distance(pts, lat, lon, ocean):
    jj, ii = np.nonzero(~ocean)
    la, lo = np.radians(lat[jj]), np.radians(lon[ii])
    for p in pts:
        a, b = np.radians(p["lat"]), np.radians(p["lon"])
        c = np.sin(a)*np.sin(la) + np.cos(a)*np.cos(la)*np.cos(lo - b)
        p["dland"] = int(round(R_EARTH * np.arccos(np.clip(c.max(), -1, 1))))


def flag_close_pairs(pts, deg=1.5):
    """Opposite-sign points within ~1.5 deg of each other are often grid noise."""
    for p in pts:
        p["pair"] = 0
    for a in range(len(pts)):
        for b in range(a + 1, len(pts)):
            p, q = pts[a], pts[b]
            if p["s"] != q["s"] and abs(p["lat"] - q["lat"]) < deg and abs(wrap(p["lon"] - q["lon"])) * np.cos(np.radians(p["lat"])) < deg:
                p["pair"] = q["pair"] = 1


def encode_field(lat, lon, amp, pha, ocean, step):
    """Complex tide on a (possibly decimated) grid, re-ordered to lon -180..180, as int16 mm (base64)."""
    a, g, oc = amp[::step, ::step], pha[::step, ::step], ocean[::step, ::step]
    lo = lon[::step]
    order = np.argsort(((lo + 180) % 360) - 180)
    a, g, oc = a[:, order], g[:, order], oc[:, order]
    zr = np.where(oc, a * 10 * np.cos(np.radians(g)), -32768)
    zi = np.where(oc, a * 10 * np.sin(np.radians(g)), -32768)
    data = np.stack([zr, zi], axis=-1).round().clip(-32768, 32767).astype("<i2")
    return {"nlat": int(a.shape[0]), "nlon": int(a.shape[1]),
            "b64": base64.b64encode(data.tobytes()).decode()}


def main(grid_dir, out_js):
    out = {"source": "GOT4.10c ocean tide, R. Ray, NASA GSFC (0.5 deg grid)", "constituents": {}}
    for c in CONSTITUENTS:
        lat, lon, amp, pha, fill = read_got_ascii(f"{grid_dir}/{c}.d")
        pts, ocean = find_points(lat, lon, amp, pha, fill)
        add_land_distance(pts, lat, lon, ocean)
        flag_close_pairs(pts)
        pts.sort(key=lambda p: -p["dland"])
        ccw = sum(p["s"] > 0 for p in pts)
        print(f"{c.upper()}: {len(pts)} points  (CCW {ccw}, CW {len(pts)-ccw}, net winding {sum(p['s'] for p in pts)})")
        for d in (0, 100, 300, 500):
            sel = [p for p in pts if p["dland"] > d and not p["pair"]]
            print(f"    > {d:3d} km from land, excluding close +/- pairs: {len(sel)}")
        print(f"    close +/- pairs (possible grid noise): {sum(p['pair'] for p in pts)//2}")
        out["constituents"][c] = {"period_h": PERIOD_H[c], "points": pts,
                                  "field": encode_field(lat, lon, amp, pha, ocean, 1)}
    with open(out_js, "w") as f:
        f.write("// Generated by find_amphidromes.py from GOT4.10c (R. Ray, NASA GSFC). Do not edit.\n")
        f.write("window.TIDE = " + json.dumps(out, separators=(",", ":")) + ";\n")
    print("wrote", out_js)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
