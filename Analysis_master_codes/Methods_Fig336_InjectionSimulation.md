# Methods

Draft methods text for the two detection-limit analyses in
`Analysis_20260506_SpatialTemporalCorrelation_Rvs.m`:

1. the analytic shot-noise detection floor (figure 336), and
2. the synthetic-event injection simulation (figures 17, 18, 22, 21).

Numbers in brackets are the parameter values as currently set in the script.
Notes for the author — inconsistencies and things to double-check before
submission — are collected at the end and are **not** part of the manuscript text.

---

## Estimation of the shot-noise-limited detection floor

### Movie preprocessing

Analysis was performed on a continuous segment of [10,001] frames
(frames [65,000–75,000]) of a motion-corrected, shutter-registered voltage
imaging movie acquired at [1 kHz] with a pixel size of [1.17 µm].

Slow photobleaching was removed before any noise estimate. The
spatially averaged fluorescence of the field of view was computed, frames
within [−5 to +10] of any detected spike were excluded, and a single
exponential was fit to the remaining time points. The fit was normalized to
its own mean and divided out of the movie frame by frame, so that the
correction is multiplicative and preserves fractional signals.

A residual movie was then formed by subtracting the temporal mean of each
pixel. Residual motion artifact was removed by projecting out, independently
at every pixel, the motion-correction traces, their squares, and the product of
the first and last motion components (linear least-squares regression, with a
constant term). All subsequent variance estimates were computed on this
motion-corrected residual movie, and the variance attributable to motion was
recorded as the difference in residual variance before and after this
regression step.

### Separating white noise from correlated signal

For each pixel we computed the total variance of the residual,

$$\mathrm{Var}_{\text{tot}}(x,y) = \big\langle\, r(x,y,t)^2 \,\big\rangle_t,$$

where $r$ is the zero-mean residual. Photon shot noise (together with read and
dark noise) is white in both space and time, whereas residual biological and
instrumental signals are correlated across neighboring pixels and across
adjacent frames. The correlated component was therefore estimated as a
covariance and subtracted, giving the white-noise (shot-noise) variance per
pixel per frame.

Two estimators were implemented. The spatial estimator, used for the analysis
shown, compares each pixel to the mean of its eight nearest neighbors,
$n(x,y,t)$, computed with a 3 × 3 kernel with the center excluded and
normalized by the number of valid neighbors so that edge pixels are treated
correctly:

$$\sigma^2_{\text{shot}}(x,y) = \mathrm{Var}_{\text{tot}}(x,y) - \big\langle\, r(x,y,t)\, n(x,y,t) \,\big\rangle_t .$$

Because shot noise is statistically independent between neighboring pixels, the
cross term $\langle r\,n\rangle$ retains only the spatially correlated
(signal) variance; subtracting it leaves the white component. The temporal
estimator is the exact analogue using the lag-one autocovariance,
$\langle r(t)\,r(t{+}1)\rangle_t$, in place of the neighbor covariance. The two
estimators are interchangeable in the code and gave consistent results.
Negative values, which arise where the estimate is noise-dominated, were
clipped to zero.

### Analytic detection floor

Analysis was restricted to pixels inside the dendritic region of interest
(the union of the first [7] dendritic footprints) with finite, positive
baseline fluorescence and a finite, positive noise variance. Baseline
fluorescence $F_0(x,y)$ was the time-averaged fluorescence after subtraction
of the camera offset ([100] counts).

Consider a hypothetical event of fractional amplitude $\alpha = \Delta F/F$
that is spatially uniform over $M$ pixels and temporally uniform over $K$
frames, detected by an unweighted (boxcar) average over those pixels and
frames. The accumulated signal is $\alpha K \sum_{i=1}^{M} F_{0,i}$, while the
noise, being white in space and time, accumulates in variance:
$K \sum_{i=1}^{M} \sigma^2_{\text{shot},i}$. Requiring a detection criterion of
$d' = $ [5] noise standard deviations gives the minimum detectable amplitude

$$\alpha_{\min}(M,K) \;=\; d'\,\frac{\sqrt{K \sum_{i=1}^{M} \sigma^2_{\text{shot},i}}}{K \sum_{i=1}^{M} F_{0,i}} \;=\; \frac{d'}{\sqrt{K}}\cdot\frac{\sqrt{\sum_{i=1}^{M} \sigma^2_{\text{shot},i}}}{\sum_{i=1}^{M} F_{0,i}},$$

which makes the $1/\sqrt{K}$ improvement with temporal averaging explicit and
reduces to $1/\sqrt{M}$ scaling in the limit of identical pixels.

Because pixels differ in brightness and noise, the order in which they are
added matters. Pixels were ranked by their individual signal-to-noise ratio,
$F_{0,i}/\sigma_{\text{shot},i}$, in descending order, so that the curve at a
given $M$ represents the best $M$-pixel unweighted average available in the
region of interest. Cumulative sums of $F_0$ and of $\sigma^2_{\text{shot}}$
along this ordering were then evaluated for $M$ up to [5000] pixels and $K$ up
to [100] frames.

For reference, the single-pixel, single-frame limit of this expression
($M = K = 1$) is reported directly, both as an absolute variance in camera
counts² and as the variance of $\Delta F/F$ for one pixel in one frame,
$\sigma^2_{\text{shot}}/F_0^2$.

Figure 336 shows $\alpha_{\min}$ on a logarithmic color scale over the range
[0.01–10] % $\Delta F/F$, with iso-detectability contours at
[0.05, 0.2 and 0.5] % $\Delta F/F$.

### Electrotonic-spread bound

Averaging over an area larger than the event itself adds noise without adding
signal, so the useful region of the (area, time) plane is bounded by how far a
real subthreshold event spreads electrotonically within the averaging window.
For a passive cable, voltage spread is formally diffusive with an effective
diffusion coefficient $D = \lambda^2/\tau_m$, where $\lambda$ is the
electrotonic length constant [110 µm] and $\tau_m$ the membrane time constant
[15 ms], giving $D \approx$ [807] µm² ms⁻¹. The axial spread after time $t$ is
$d = \sqrt{2Dt}$, and for bidirectional spread along a dendrite of width
[5 µm] the corresponding membrane area is $A = 2\,d\,w$, converted to a pixel
count using the measured pixel size. This bound is overlaid on figure 336
(dashed line): points to the right of the line correspond to spatial averaging
windows larger than the event itself at that temporal averaging window.

---

## Detection limit measured by synthetic-event injection

The analytic floor above assumes an ideal, matched observer. To measure the
detection limit actually achieved by our analysis pipeline, we injected
synthetic subthreshold events of known amplitude, spatial extent and time
course into real imaging data and recovered them with the unmodified detection
code.

### Background data

Injections were performed on the same motion-corrected residual movie, so that
each trial contains the true shot noise, the true residual motion artifact and
the true spontaneous activity of the recording. For each condition, [40]
background windows of [1001] frames were drawn at random (with a fixed random
seed for reproducibility) from the [10,001]-frame segment.

### Synthetic events

A synthetic event was added to the residual movie as a separable
space–time perturbation of the local baseline fluorescence,

$$I_{\text{inj}}(x,y,t) \;=\; I_{\text{base}}(x,y,t) \;+\; s\,\alpha\, F_0(x,y)\, S(x,y)\, h(t),$$

with $s = -1$ to match the negative-going fluorescence response of Voltron2.
Scaling by the local $F_0$ makes $\alpha$ a genuine fractional change
($\Delta F/F$) at every pixel rather than a fixed number of counts.

**Spatial profile.** $S(x,y)$ was constructed on the traced dendritic
skeleton. Skeleton nodes were weighted by a Gaussian in arc length from an
anchor node, with standard deviation $\sigma_L = \ell/2$, and each weighted
node was rendered as a two-dimensional Gaussian of width [2 px] ([2.3 µm]);
the result was normalized to unit peak. The injected event therefore follows
the dendrite rather than being an isotropic blob, and $\sigma_L$ (the
half-extent, plotted on the abscissa of figures 17, 18 and 22) parameterizes
its spatial extent. Full extents of $\ell = $ [1 to 381 µm in 20 µm steps]
were tested.

**Time course.** $h(t)$ was a peak-normalized alpha function,
$h(t) = (t/\tau)\,e^{\,1-t/\tau}$ for $t > 0$, with onset at the center of the
window and peak at $t_{\text{onset}} + \tau$. Time constants of $\tau = $
[1 to 66 ms in 5 ms steps] were tested.

**Amplitude.** $\alpha = $ [0.5, 1, 2, 4, 8 and 14] % $\Delta F/F$. For
comparison with electrophysiology, amplitudes were also expressed in millivolts
using the measured sensitivity of [12 % $\Delta F/F$ per 100 mV].

Events were injected at [3] anchor positions along the dendrite, one drawn at
random on each repetition, so that the reported probabilities average over
local variation in brightness and noise. In total, [20] spatial extents ×
[6] amplitudes × [14] time constants = [1680] conditions were each tested
[40] times.

### Event recovery

Injected movies were analyzed exactly as experimental data. Fluorescence
traces were extracted by footprint-weighted projection of the movie onto each
dendritic region of interest, sign-inverted, and normalized by a slowly varying
fit to the somatic spike amplitude so that traces are expressed in spike-height
units and are insensitive to residual bleaching. Detection was scored
separately for spikes and for subthreshold peaks.

**Spikes.** Each trace was divided by a robust noise estimate obtained from the
median-filtered trace using only its negative-going half (mirrored to form a
symmetric distribution), making the estimate insensitive to positive-going
events. Spikes were then detected with a threshold of [2.5] and a local
prominence of [1.5] noise standard deviations, the prominence being evaluated
on a [7]-point high-pass-filtered version of the trace. To count only
injection-evoked spikes, spikes belonging to the recorded spike train and
spikes independently detected in the same background window *before* injection
were both subtracted. An injection was scored as detected if a residual spike
fell within [±7] frames of the injected waveform peak.

**Subthreshold peaks.** Traces were converted to $\Delta F/F$ using the
per-region scale factor calibrated from the mean back-propagating action
potential. Spikes were masked with a [±3]-frame dilation, the trace was
linearly interpolated across the masked samples and smoothed over [20] frames
to isolate the subthreshold component, and the result was normalized by the
mean back-propagating action potential amplitude measured from the
spike-triggered average of the same neuron. Peaks were detected with a minimum
height of [0.2] and a minimum prominence of [0.1] spike-height units and a
minimum separation of [50] frames. An injection was scored as detected if a
peak fell within the support of the injected waveform (defined as
$h(t) > 0.1$) in a region of interest overlapping the injection site.

Regions of interest were classified per anchor: the *injected* set comprised
footprints containing the anchor pixel, and the *checked* set additionally
included all footprints within [35 µm] measured along the dendritic tree.
Scoring peaks and spikes only within the injected set makes the reported
probabilities spatially specific.

### Reported quantities

For each condition, the detection probability was the fraction of the [40]
repetitions in which the event was recovered, computed separately for
subthreshold peaks, for spikes, and for their per-repetition union
(peak **or** spike). Because the union is evaluated per repetition rather than
from the two marginal probabilities, it is the true joint detection
probability.

Figure 17 shows peak and spike detection probability as a function of spatial
extent $\sigma_L$ and time constant $\tau$, one column per amplitude, with the
electrotonic-spread bound overlaid. Figure 18 shows the same data as a
function of amplitude and spatial extent, one column per time constant.
Figure 22 shows the probability of detecting the event by either criterion, one
map per amplitude, on the axes of figure 17. Figure 21 shows detection
probability against $\tau$ at a fixed spatial extent of [56 µm], for each
amplitude.

---

## Notes for the author (not manuscript text)

1. **Camera offset is subtracted inconsistently between the two analyses.**
   The noise-floor section uses $F_0 = \langle I\rangle_t - 100$
   (`Analysis_20260506_SpatialTemporalCorrelation_Rvs.m:429`), but the
   injection section uses $F_0 = \langle I\rangle_t$ with the offset left in
   (line 622). Injected events are therefore scaled by $F + 100$ rather than
   $F$, making the true injected $\Delta F/F$ slightly *smaller* than the
   nominal $\alpha$ — the measured detection thresholds are correspondingly
   conservative. Worth fixing, or at least quantifying, since the two figures
   are presented as complementary.
2. **Abscissa units in figure 336.** The x axis is labeled
   "Spatial averaging area (µm²)" but the plotted variable is a pixel count.
   At 1.17 µm per pixel, one pixel is 1.37 µm², so the axis is off by that
   factor if read as an area. Either relabel to "pixels averaged" or multiply
   `M_list` by `pix_um^2`.
3. **"Shot noise" is really "white noise."** The estimator captures all
   spatially and temporally uncorrelated variance, which includes read and dark
   noise. The methods text above says this; if you prefer the shorter term in
   the figure and main text, consider defining it once.
4. **Uniform weighting is not optimal.** Ranking pixels by SNR and then
   averaging them with equal weights is a defensible practical bound, but a
   matched filter (weights $\propto F_0/\sigma^2$) would do better. Worth one
   sentence stating that $\alpha_{\min}$ is the floor for uniform weighting,
   not the Cramér–Rao-style optimum.
5. **No explicit false-positive rate for the injection analysis.** The
   analytic floor uses $d' = 5$, whereas the injection analysis reports hit
   rates from fixed detector thresholds. Running the detector on
   non-injected windows to quote a false-positive rate per window would make
   the two panels directly comparable.
6. **Syntax error at line 419.** A stray `ƒ` character precedes the cell
   marker (`ƒ%% Calculate low bound threshold ∆F/F`); the file will not run
   top-to-bottom as a script until it is removed.
