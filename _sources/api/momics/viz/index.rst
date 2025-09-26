momics.viz
==========

.. py:module:: momics.viz


Attributes
----------

.. autoapisummary::

   momics.viz.DEFAULT_ISM_CMAP


Functions
---------

.. autoapisummary::

   momics.viz.aggrcoverage
   momics.viz.coverage
   momics.viz.heatcoverage
   momics.viz.plot_ISM_heatmap


Module Contents
---------------

.. py:function:: aggrcoverage(q, ci = True, cmap = None, fig_size = None)

   Plot the aggregate coverage of tracks queried over a number of genomic windows.

   :param q: A MomicsQuery object containing the coverage data.
   :type q: MomicsQuery
   :param ci: Whether to plot confidence intervals. Defaults to True.
   :type ci: bool, optional
   :param cmap: A dictionary mapping track names to colors. Defaults to None.
   :type cmap: dict, optional
   :param fig_size: Size of the figure. Defaults to (10, 1.5 * n_tracks).
   :type fig_size: Tuple[int, int], optional

   :returns: The matplotlib Axes object containing the plot.
   :rtype: plt.Axes


.. py:function:: coverage(q, cmap = None, fig_size = None)

   Plot the coverage of tracks queried over a single genomic window.

   :param q: A MomicsQuery object containing the coverage data.
   :type q: MomicsQuery
   :param cmap: A dictionary mapping track names to colors. Defaults to None.
   :type cmap: dict, optional
   :param fig_size: Size of the figure. Defaults to (10, 1.5 * n_tracks).
   :type fig_size: Tuple[int, int], optional

   :returns: The matplotlib Axes object containing the plot.
   :rtype: plt.Axes

   .. rubric:: Example

   >>> repo = Momics("tests_data/S288c_MTL.momics")
   >>> q = MomicsQuery(repo, "IV:324469-332598").query_tracks(tracks=["MNASE", "H3K4ME3", "ATAC", "SUA7"])
   >>> fig = momics.viz.coverage(q)


.. py:function:: heatcoverage(q, order_by = None, cmap = None, fig_size = None)

   Plot heatmaps of queried coverage data over genomic regions.
   Each track is plotted in a separate heatmap.

   :param q: A MomicsQuery object containing the coverage data.
   :type q: MomicsQuery
   :param order_by: Track name to order the regions by. Defaults to None.
   :type order_by: str, optional
   :param cmap: A dictionary mapping track names to colors. Defaults to None.
   :type cmap: dict, optional
   :param fig_size: Size of the figure. Defaults to (10, 1.5 * n_tracks).
   :type fig_size: Tuple[int, int], optional

   :returns: The matplotlib Axes object containing the plot.
   :rtype: plt.Axes


.. py:function:: plot_ISM_heatmap(ism_pos, figsize=(50, 2), cmap=DEFAULT_ISM_CMAP, figname='tmp.pdf')

   Plot a heatmap of mutation impacts.
   :param ism_pos: Dictionary with position as key and nucleotide impacts as values
   :param figsize: Size of the figure
   :param cmap: Colormap for the heatmap
   :param figname: Filename to save the figure


.. py:data:: DEFAULT_ISM_CMAP
   :value: 'afmhot_r'


