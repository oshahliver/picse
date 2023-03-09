#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 30 23:15:52 2019

@author: oshah
"""

# from picse.materials import brucite_phase as bruce
from picse.materials import phase_transitions_water_Wagner2002 as water
from picse.physicalparams import phase_transition_coeffs_MgSiO3
from picse.materials import eos


def getPhase(ll=0, T=None, P=None, d=None, double_check=False):

    if d == None and double_check:
        d = eos.Compute(what="dens", T=T, P=P, ll=ll)

    # Water
    if ll == 0:
        return water.phaseCheck(P, T, double_check=False)

    # Magnesiowüstite
    elif ll == 11:
        raise NotImplementedError("Mg(OH)2 EoS is currently not available here.")
        return bruce.phaseCheck(T=T, P=P)

    # Perovskite
    elif ll == 6:
        P0, a = phase_transition_coeffs_MgSiO3

        # Compute transition pressure in GPa
        Ptrans = P0 + a * T

        # Convert to Pa
        Ptrans *= 1.0e9

        if P < Ptrans:
            return 0

        else:
            return 1
