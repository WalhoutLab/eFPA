# -*- coding: utf-8 -*-
import cobra

# Load model
model = cobra.io.read_sbml_model('../ModelFiles/xml/yeastGEM.xml')

# Correct metabolite ids:
for met in model.metabolites:
    met.id = met.id.replace('__91__', '_')
    met.id = met.id.replace('__93__', '')

