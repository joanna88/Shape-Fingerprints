#!/usr/bin/env python3

from __future__ import print_function

from openeye.oechem import *
from openeye.oeshape import *
import time

def comp(refmol, fitmol):

    """
    additional function in case the fitting molecule has bigger volume
    - enable to compare two molecules instead of whole set
    """

    best = OEBestOverlay()
    best.SetRefMol(refmol)

    #best.SetInitialOrientation(OEBOOrientation_Subrocs)
    best.SetInitialOrientation(OEBOOrientation_InertialAtHeavyAtoms)

    scoreiter = OEBestOverlayScoreIter()
    OESortOverlayScores(scoreiter, best.Overlay(fitmol), OEHighestTanimoto())

    return scoreiter

def get_fingerprint(mol, Shape_database, bitOn):
    """
    Calculates Shape Tanimoto for molecule in reference to molecules stored in Shape Database,
    returns on bits for value higher than Design Tanimoto
    and off bit for value smaller
    :param DT: bitOn value (threshold value)
    :return: fingerprints for input molecule
    """

    keepsize = 1

    best = OEBestOverlay()
    best.SetRefMol(mol)

    # sets the Initial orientation of molecules - subrocs

    best.SetInitialOrientation(OEBOOrientation_InertialAtHeavyAtoms)
    #best.SetInitialOrientation(OEBOOrientation_Subrocs)

    V_ref = OECalcVolume(mol)

    fingerprint = ''

    for fitmol in Shape_database.GetOEMols():

        #print(fitmol.GetTitle())

        resCount = 0

        V_fit = OECalcVolume(fitmol)

        # molecules are all the time compared the same way
        # - molecule with bigger volume is always the reference

        if V_fit > V_ref:

            scoreiter = comp(fitmol, mol)

        else:

            scoreiter = OEBestOverlayScoreIter()
            OESortOverlayScores(scoreiter, best.Overlay(fitmol), OEHighestTanimoto())

        for score in scoreiter:
            outmol = OEGraphMol(fitmol.GetConf(OEHasConfIdx(score.fitconfidx)))
            score.Transform(outmol)
            #print(score.tanimoto)

            if float(score.tanimoto) > bitOn:
                fingerprint += ' 1'
            else:
                fingerprint += ' 0'

            resCount += 1
            if resCount == keepsize:
                break

    return fingerprint

def main(argv=[__name__]):

    if len(argv) != 5:
        OEThrow.Usage("%s <data_file.sdf> <shape_file.sdf> <out_file.sdf> <bitOn> " % argv[0])

    data_file = oemolistream(argv[1])
    shape_file = OEMolDatabase(argv[2])
    out_file = oemolostream(argv[3])
    bitOn = float(argv[4])

    start = time.time()
    mol_count = 0

    for mol in data_file.GetOEMols():
        fp = get_fingerprint(mol, shape_file, bitOn)
        OESetSDData(mol, "Fingerprint", "%s" % fp)
        OEWriteMolecule(out_file, mol)
        mol_count += 1

    end = time.time()
    print("total time: ", end-start, "   time per molecule:", (end - start)/mol_count)


if __name__ == "__main__":
    import sys
    sys.exit(main(sys.argv))
