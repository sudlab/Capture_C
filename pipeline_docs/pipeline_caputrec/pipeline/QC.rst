Quality Control
================


Read categoriation
-------------------

.. report:: QC.CountMetrics
   :render: r-ggplot
   :transform: melt
   :groupby: all
   :statement: aes(x=track, y=value/1000000, fill=variable) + geom_bar(position="stack", stat="identity") + xlab("Track") + ylab("Reads (millions)") + theme_bw() + scale_fill_discrete(name="Category")

   Number of reads for each library with correct mapping geometry 


Library anomolies
--------------------


.. report:: QC.Anomolies
   :render: r-ggplot
   :groupby: all
   :statement: aes(x=slice, y=count, fill=category) + geom_bar(position="stack", stat="identity") + theme_bw() + theme(axis.text.x=element_text(angle=-90)) + xlab("Probe") + ylab("Reads") + facet_grid(track~.)

   Number of reads mapping to or adjecent to probe fragment

