from CGATReport.Tracker import *
from CGATReport.Utils import PARAMS as P

PIPELINE_DB = P.get("capturec_database",
                    P.get("database_name",
                          None))

class CaptureCTracker(TrackerSQL):
    
    def __init__(self, *args, **kwags):

   	print "new"
	print PIPELINE_DB
        TrackerSQL.__init__(self, backend=PIPELINE_DB,
                            *args,
                            **kwargs)

    
