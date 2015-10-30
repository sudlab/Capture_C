from CaptureCTracker import CaptureCTracker
from CGATReport.Tracker import SQLStatementTracker


class CountMetrics(CaptureCTracker, SQLStatementTracker):

    statement = '''SELECT track, Unpaired, no_probe, multi_hit,
                  fragments - (Unpaired + no_probe + multi_hit) as OK
                   FROM interaction_count_metrics'''
    

class Anomolies(CaptureCTracker):

    def getTracks(self):

        return self.getValues("SELECT DISTINCT track FROM interaction_counts")

    def getSlices(self):

        return self.getValues("SELECT DISTINCT probe FROM probe_fragments_lookup")


    def __call__(self, track, slice=None):

        statement = '''SELECT DISTINCT track, probe, count, frag1, frag2
                       FROM interaction_counts as ic
                         INNER JOIN probe_fragments_lookup as lu
                          ON ic.frag1 = lu.fragment
                       WHERE track = '%(track)s' AND probe = '%(slice)s'
                    '''

        result = self.getDataFrame(statement % locals())

        self_lig =  result["Count"][result["Frag1"] == result["Frag2"]].sum()
        undigested =  result["Count"][abs(result["Frag1"] - result["Frag2"]) == 1].sum()
        other = result["Count"].sum() - self_lig - undigested

        return {"category": ["Self ligation", "Undigested", "other"],
                "count": (self_lig, undigested, other)}
