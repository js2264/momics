import os
import tiledb
import numpy as np
import pandas as pd

__all__ = ["Momics"]

class Momics:
    """
    A convenient interface to a momics data collection.

    Parameters
    ----------
    store : str, Path to a momics file (URI string)
    kwargs : Options passed to methods.

    Notes
    -----
    If ``store`` is a file path, the file will be opened temporarily in
    when performing operations. This allows :py:class:`Momics` objects to be
    serialized for multiprocess and distributed computations.

    Table selectors, created using :py:meth:`genome`, :py:meth:`scores`, and
    :py:meth:`intervals`, perform range queries over table rows,
    returning :py:class:`pd.DataFrame`, :py:class:`pd.Series` or
    :py:class:`numpy.ndarray`

    Metadata is accessible as a dictionary through the :py:attr:`info`
    property.

    """

    def __init__(self, store: str):
        self._refresh()

    def _refresh(self) -> None:
        try:
            None
        except KeyError:
            raise KeyError("err_msg") from None

    def genome(self, **kwargs):
        """Genome array selector

        Returns
        -------
        Array

        """

        return self

    def scores(self, **kwargs):
        """Scores array selector

        Returns
        -------
        Array

        """

        return self

    def intervals(self, **kwargs):
        """Intervals array selector

        Returns
        -------
        Array

        """

        return self

    # @property
    # def info(self) -> dict:
    #     """File information and metadata
    #
    #     Returns
    #     -------
    #     dict
    #
    #     """
    #     with open_hdf5(self.store, **self.open_kws) as h5:
    #         grp = h5[self.root]
    #         return info(grp)

