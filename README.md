# geodetic_inversion
This repo homogeneous/layered inversion using InSAR/GPS

---
## Step 1 ~ 2 are written in the file main_detrend_inversion.m
### Step 0: setup your MATLAB, CSHELL and GMT paths.

### Step 1: data cleaning using clean_insar_data.m
```MATLAB
% remove some near-field unwrapping errors manually first
clean_insar_data; 
```
If you open `clean_insar_data.m`, you need to specify following information to clean the data:
- `this_track`: the directory that InSAR grid file is saved.
- `insar_file = 'unwrap_ll.grd'`: the interferogram/offsets grid
- `mask_file = 'mask_txt'`: the masking TXT file that writes the masking polygon files (e.g. mask1.txt, mask2.txt, ...). 
The masking polygons could be generated using Google Earth KML files.

Then the subroutine runs as (inside clean_insar_data.m):
``` MATLAB
% scale = -lambda / 4 / pi
% scale = -4*pi if the insar_file is an offset grid
% los_max is used just for plot (cm)
% detrend = 1 means that you apply for a detrend with the topography
mask_insar_phase(this_track, insar_file, mask_file, scale, 'los_max', 80, 'detrend', 0);
```
This step would generate a subsampled grid file called "unwrap_clean_sample.grd", in 100 meters resolution as default.

### Step 2: detrend the phase and remove the phase ambiguity
If we have enough far-field GPS data, we could use those GPS data to invert a coarse slip model to detrend the unwrapped phase. \
In cases such as Pamir and Qinghai earthquake, since we do not have enough GPS sites covered, we could just assume a far-field pixel that corresponds to zero displacement.
```MATLAB
% grdin = '/Users/zej011/coseismic/DES5/LOS3/unwrap_clean_sample.grd';
% grdout = '/Users/zej011/coseismic/DES5/LOS3/los_clean_detrend.grd';
% lonf = 73.215150;   latf = 37.754558;  ref_lon = 71;   
% threshold = 0.5; (pixels within 500m would be averaged in order to get the value at the reference point)
remove_ref_from_grid(grdin,grdout,lonf,latf,ref_lon,threshold);
```
- `grdin`: InSAR grid input
- `grdout`: InSAR grid output
- `lonf/latf`: reference lon/lat to convert lon/lat to UTM coordinates (m or km)
- `ref_lon`: central meridian to compute UTM coordinates: This parameter is necessary because UTM coordinates are generally confined within a single UTM zone.
If your study area crosses two different UTM zones, it would be complicated to convert to the UTM coordinates at the UTM zone boundary directly. 
So we shift the central meridian to the west (ref_lon), so that you could define a broader UTM zone than the common one. Usualy choose a ref_lon <= lonf.

### (Optional) apply the sign mask for the detrended offset data
Because offset data are noisy than phase ones, sometimes you even need to apply a sign mask across the fault.
In each directory (this_track), just put two polygons named with clean_left.txt and clean_right.txt that surround the left and right part of the offset map.
```MATLAB
cd(this_track);
movefile los_clean_detrend.grd los_clean_unmask.grd
sign_mask_offset(this_track, 'los_clean_unmasked.grd');
```
