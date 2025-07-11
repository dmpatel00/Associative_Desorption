import taskblaster as tb
from pathlib import Path
import runpy

tasks_path = Path(__file__).parent / 'tasks2.py'
globals().update(runpy.run_path(tasks_path))

@tb.workflow
class ParametrizedAdsorbateWorkflow:
    calculator = tb.var()

    @tb.dynamical_workflow_generator({'symbols': '*/surf'})
    def surfaces(self):
        return tb.node('parametrize_adsorbate_workflow',
                       calculator=self.calculator)


def workflow(rn):
    calculator = {
        'mode': {'name':'pw',
            'ecut':500},
        'kpts': {'size': (4,3,1)},
        'xc': 'BEEF-vdW',
        'txt': 'gpaw.txt'
    }
    wf = ParametrizedAdsorbateWorkflow(calculator=calculator)
    rn.run_workflow(wf)
