import numpy as np
from scipy.stats import gaussian_kde, norm
from scipy.signal import find_peaks

import matplotlib as mpl; mpl.use("PDF")
import matplotlib.pyplot as plt
import seaborn as sns

#-------------------------------- Content Class --------------------------------#

class _Peaks:

    def __init__(self, x, y, peaks, props):
        self.values = x        # original data of distribution
        peaks = peaks[np.nonzero((0.02 < x[peaks]) & (x[peaks] < 0.48))]
        self.index = peaks     # index of peaks
        self.x = x[self.index] # X coordinate of peaks (MAF)
        self.y = y[self.index] # Y coordinate of peaks (Density)
        self.left_ips = np.round(props['left_ips']).astype(int)
        self.right_ips = np.round(props['right_ips']).astype(int)
        # Interpolated positions of left and right intersection points of a 
        # horizontal line at the respective evaluation height.

    def get_density_matrix(self):
        """
        MAF distribution is composed of several local MAF distributions,
        which are assumed to be normal distributions with `$mu$ = peak` and
        `$sigma$ = abs(left_ips - right_ips).`
        """
        data = self.values
        l_ips = self.values[self.left_ips]
        r_ips = self.values[self.right_ips]

        density_matrix_mix = []
        # shape: n_peaks (rows) x n_peaks (columns)
        # each row is an array of estimated density of the peak in another
        # location of peaks.
        for mu, scale, l, r in zip(self.x, self.y, l_ips, r_ips):
            _, sigma = self._single_peak_normal_distribution(data, l, r)
            density_matrix_mix.append(
                self._single_peak_scaled_density(data, mu, sigma, scale))
        else:
            density_matrix_mix = np.round(np.array(density_matrix_mix), 4)
            # density_matrix.shape == (len(peaks), len(peaks))

        # Left peaks: Fetal specific allele AAab
        # Right peaks: Maternal specific allele ABaa
        AAabDM = density_matrix_mix[np.where(self.x <= .25)].sum(axis=0)
        ABaaDM = density_matrix_mix[np.where(self.x > 0.25)].sum(axis=0)
        AAabDM = self.y * (AAabDM / (AAabDM + ABaaDM))
        ABaaDM = self.y * (ABaaDM / (AAabDM + ABaaDM))
        return AAabDM, ABaaDM

    def _single_peak_normal_distribution(self, data, left_bound, right_bound):
        index = np.where((left_bound <= data) & (data <= right_bound))
        return norm.fit(data[index]) # mu and sigma of a single peak

    def _single_peak_scaled_density(self, data, mu, sigma, scale):
        den = norm.pdf(data, loc=mu, scale=sigma) # probability density function
        return den[self.index] / den.max() * scale

class MafDistribution:

    def __init__(self, depths=None, values=None, *,
                 bw_method=0.05, depth_weighted=None):
        self.depths = depths or []
        self.values = values or []

        if depth_weighted:
            depth_weighted = depth_weighted if callable(depth_weighted) \
                        else lambda arr: np.log1p(arr) + 1
            self.depths = depth_weighted(np.array(self.depths)).round()
            self.values = np.repeat(self.values, self.depths.astype(int))

        self.kernel = self._gaussian_kde(bw_method=bw_method)

        if self.kernel is None:
            self.nloci = 0
        else:
            self.nloci = self.kernel.neff
            self.peaks = None
            self.distribution = None
            # self.peaks and self.distribution will be filled after
            # `find_peaks` method is called.

    def depth_sum(self):
        return sum(self.depths)

    def _gaussian_kde(self, *, bw_method):
        data = np.array(self.values)
        data = data[np.where((0 <= data) & (data <= 0.5))]
        try:
            return gaussian_kde(data, bw_method=bw_method)
        except np.linalg.LinAlgError:
            return None

    def find_peaks(self, *, width=0, rel_height=1, **kwargs):
        x = np.linspace(0, 0.501, 502*10-9) # accuracy: 1e-4
        y = self.kernel.evaluate(x)
        y = y / y.max() # normalize to 0 ~ 1
        peaks, props = find_peaks(y, width=width, rel_height=rel_height, **kwargs)
        setattr(self, 'distribution', (x, y))
        setattr(self, 'peaks', _Peaks(x, y, peaks, props))

    def estimate_fraction(self):
        AAabDM, ABaaDM = self.peaks.get_density_matrix() # DM: density matrix
        peak_AAab = self.peaks.x.dot(AAabDM) / AAabDM.sum()
        peak_ABaa = self.peaks.x.dot(ABaaDM) / ABaaDM.sum()
        frac_AAab = 2 * peak_AAab
        frac_ABaa = 1 - 2 * peak_ABaa
        frac_wAvg = (frac_AAab*AAabDM.sum() + frac_ABaa*ABaaDM.sum()) \
                 / (AAabDM.sum() + ABaaDM.sum())

        fractions = np.array([frac_AAab, frac_ABaa, frac_wAvg])

        if np.isnan(fractions).any():
            not_nan_frac = fractions[np.where(~np.isnan(fractions))]
            not_nan_frac = np.squeeze(not_nan_frac)
            frac_AAab = frac_ABaa = frac_wAvg = not_nan_frac

        return frac_AAab, frac_ABaa, frac_wAvg

    def graph(self, outfig=None):
        BLACK, BLUE, RED = 'k', '#1867bb', '#d3121b'
        x, y = self.distribution
        peak_x = self.peaks.x
        peak_y = self.peaks.y

        fig, ax = plt.subplots(figsize = (12, 8))

        #ax.hist(self.values, bins=50, color='k', alpha=.3, density=True)
        ax.plot(x, y, color=BLACK, linewidth=2,
                label='MAF distribution (smoothed by Gaussian kernel)')

        for h, (i, j) in enumerate(zip(peak_x, peak_y)):
            color = BLACK if h in (0, len(peak_x)-1) \
               else BLUE if i < 0.25 \
               else RED
            scl_x = ax.get_xlim()[1]
            scl_y = ax.get_ylim()[1]
            txt_i = f'maf: {np.round(i, 4)}'
            txt_j = f'den: {np.round(j, 4)}'
            ax.annotate(txt_i, xy=(i, j), xytext=(i, j+scl_y*.045), ha='center', color=color)
            ax.annotate(txt_j, xy=(i, j), xytext=(i, j+scl_y*0.02), ha='center', color=color)
            ax.axvline(i, 0, j/scl_y, color=color, linestyle='--')

        ax.set(
            xlim = (x.min(), x.max()),
            yticks = [],
            xlabel = 'minor allele frequency',
            ylabel = 'density',
        )
        ax.legend()
        fig.savefig(outfig, dpi=800, transparent=True, bbox_inches='tight')
