import taskblaster as tb


@tb.workflow
class ParametrizedMaterialsWorkflow:
    calculator = tb.var()

    @tb.dynamical_workflow_generator({'symbols': '*/atoms'})
    def systems(self):
        return tb.node('parametrize_materials_workflow',
                       calculator=self.calculator)


def workflow(rn):
    calculator = {
        'mode': {'name':'pw',
            'ecut':500},
        'kpts': {'size': (8,8,8)},
        'xc': 'BEEF-vdW',
        'txt': 'gpaw.txt'
    }
    wf = ParametrizedMaterialsWorkflow(calculator=calculator)
    rn.run_workflow(wf)
