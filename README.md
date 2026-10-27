# Amphidromic Earth

An interactive 3D globe of Earth's **amphidromic points**: the places in the ocean where a tide has zero height and high water circles around instead of rising and falling.

**Open it:** https://krimseylilleth.github.io/amphidromic-earth/

## What you can do

- Watch the real tide move around the globe for the four largest tides: M2 (main lunar, 12.42 h), S2 (solar, 12 h), K1 and O1 (once a day).
- See every amphidromic point found in the data, colored by which way the tide turns around it.
- Click anywhere on the ocean to get a side view through that spot. It shows:
  - the tide surface
  - the real sea floor
  - the water currents
  - a current ellipse (the current's path over one cycle, seen from above)
- Read live check numbers next to every picture, recomputed in the page from the data.

## Where the data comes from

| What | Source |
|---|---|
| Tide heights and timing | **GOT4.10c** ocean tide model, R. Ray, NASA Goddard Space Flight Center, built from satellite altimetry, 0.5° grid. [Download page](https://earth.gsfc.nasa.gov/geo/data/ocean-tide-models) |
| Sea floor depth | **NOAA ETOPO1**, sampled every 10 arc-minutes via [NOAA ERDDAP](https://coastwatch.pfeg.noaa.gov/erddap/griddap/etopo180.html), averaged to 0.5° |

## How the points are found

The tide at each grid point is `A cos(wt - G)`, where `A` is the amplitude and `G` the phase lag. The script `scripts/find_amphidromes.py` walks around every 0.5° grid cell and adds up the change in `G`.

- **No full turn:** if `G` doesn't wind a full 360°, there is no point in the cell.
- **Full turn:** if `G` winds a full 360°, the cell contains an amphidromic point, a spot where the amplitude is zero. Winding one way means the tide turns counter-clockwise (seen from above); winding the other way means clockwise.

The exact position inside the cell comes from finding the zero of the interpolated tide.

**Counts for M2:**

| | Points |
|---|---|
| Found in the model | 74 |
| Shown by default (tight pairs hidden) | 65 |
| At least 300 km from land | 12 |

The usual rule (counter-clockwise in the northern hemisphere, clockwise in the southern) holds for 82% of them.

## How the currents are worked out

The model gives only water height, so the currents are calculated. The slope of the tide surface pushes the water, and Earth's rotation turns it sideways. This is the linear shallow-water momentum balance, with no friction.

This works well in the deep ocean, where currents are about 1 to 3 cm/s. In shallow seas the real currents are slower because of friction, so the numbers there are upper limits. Near the latitude where a tide's period matches the turning period of Earth's rotation (74.5° for M2, about 30° for K1 and O1), the calculation breaks down and the page doesn't show a current.

## Limits

- **The grid is coarse.** A 0.5° grid (about 55 km) can't see points in narrow seas, straits and bays. It also misses points that sit on land, such as near southern Norway, eastern Ireland and New Zealand. Finer tide models find more points near coasts.
- **The animation time is arbitrary.** It starts at an arbitrary moment in the cycle, so it isn't today's actual tide. Times are correct relative to each other.
- **Each tide has its own points.** The real sea surface is the sum of many tides, and this globe shows the four largest.

## Suggest an improvement

- **Ideas, questions, errors:** open an [Issue](https://github.com/KrimseyLilleth/amphidromic-earth/issues).
- **Code changes:** fork the repository and send a pull request.

## Rebuild the data yourself

You need Python 3 with numpy and scipy installed.

```
# 1. Download and unpack GOT4.10c
curl -L -o got4.10c.tar.gz https://earth.gsfc.nasa.gov/sites/default/files/2023-12/got4.10c.tar.gz
tar -xzf got4.10c.tar.gz
python scripts/find_amphidromes.py GOT4.10c/grids_oceantide amphidromic_got410c_data.js

# 2. Download ETOPO1 at 10 arc-minute spacing
curl -o etopo_10min.nc "https://coastwatch.pfeg.noaa.gov/erddap/griddap/etopo180.nc?altitude%5B(-90.0):10:(90.0)%5D%5B(-180.0):10:(180.0)%5D"
python scripts/make_bathy_grid.py etopo_10min.nc amphidromic_bathy_data.js
```

## Files

| File | What it is |
|---|---|
| `index.html` | The globe (three.js) |
| `amphidromic_got410c_data.js` | Tide fields and points, generated from GOT4.10c |
| `amphidromic_bathy_data.js` | Sea floor depth, generated from ETOPO1 |
| `scripts/find_amphidromes.py` | Finds the points and writes the tide data file |
| `scripts/make_bathy_grid.py` | Averages ETOPO1 onto the tide grid |

---

Made by Krimsey Lilleth, klilleth@proton.me
