import numpy as np
from scipy.signal import fftconvolve, find_peaks

class SignalArray(np.ndarray):

    def __new__(cls, array_obj, smoothed=False):
        ins = np.asarray(array_obj).view(cls)
        ins.smoothed = smoothed
        return ins

    def smoothen(self, kernel_size, mode='valid'):
        # smoothen ns by average
        #self.smoothed = True
        kernel = np.ones(kernel_size) / kernel_size
        return type(self)(fftconvolve(self, kernel, mode=mode), smoothed=True)

    def normalize(self, method='minmax', mean=None):

        if method == 'minmax':
            self._minmax_normalize()

        elif method == 'mean':
            self._mean_normalize(mean=mean)

        else:
            options = ("minmax", "mean")
            raise ValueError(f'Valid options for `method`: {options}.')

    def find_peaks(self, distance=147, prominence=0.2, **kwargs):
        if not self.smoothed:
            raise SyntaxError("Method `find_peaks` should be called after array "
                              "had been smoothed by method `smoothen`.")
        kwargs.update(distance=distance, prominence=prominence)
        peaks_up,   _ = find_peaks( self, **kwargs)
        peaks_down, _ = find_peaks(-self, **kwargs)
        return peaks_up, peaks_down

    def _minmax_normalize(self):
        max_ = self.max()
        min_ = self.min()
        diff = max_ - min_
        self -= min_
        self /= diff

    def _mean_normalize(self, mean=None):
        if mean is None:
            self /= self.mean()
        else:
            self /= mean
