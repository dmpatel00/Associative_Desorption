import taskblaster as tb

@tb.workflow
class ParametrizedEndStateWorkflow:
    calculator = tb.var()

    @tb.dynamical_workflow_generator({'symbols': '*/surf'})
    def endstates(self):
        return tb.node('parametrize_endstate_workflow',
                       calculator=self.calculator)


def workflow(rn):
    calculator = {
        'mode': {'name':'pw',
            'ecut':500},
        'kpts': {'size': (4,3,1)},
        'xc': 'BEEF-vdW',
        'txt': 'gpaw.txt'
    }
    wf = ParametrizedEndStateWorkflow(calculator=calculator)
    rn.run_workflow(wf)
