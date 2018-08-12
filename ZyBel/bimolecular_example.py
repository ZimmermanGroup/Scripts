import zstruct2
import pybel
import argparse
import confab

# => Please see license.txt for licensing and copyright information <= 
#       =>    Zimmerman Group, University of Michigan <= #


def main(substrate1,substrate2,read_xyz=False):

    smiles1=substrate1
    smiles2=substrate2

    # ============== Generate Molecule ========== #
    if read_xyz:
        mol1=pybel.readfile("xyz","react1.xyz").next()
        mol2=pybel.readfile("xyz","react2.xyz").next()
    else:
        mol1=confab.gen3d(smiles1)
        mol2=confab.gen3d(smiles2)

    # ============ Driving Coordinate Idx Generation ==============#
    # =======> Smarts Patterns used to find reactive indices <=====#
    smarts1 = pybel.Smarts(substrate1)
    smarts2 = pybel.Smarts(substrate2)
    mol1_idx=smarts1.findall(mol1)
    mol2_idx=smarts2.findall(mol2)
    r1=[]
    r2=[]
    for i in mol1_idx:
        r1+=list(i)
    for i in mol2_idx:
        r2+=list(i)
    for i in range(len(r1)):
        r1[i]+=len(mol1.atoms)

    #==================== do zstruct ============================ #
    zstruct2.zstruct2(mol1,mol2,r1,r2,None,doOne=True,doTwo=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parse Bool")   
    parser.add_argument('--substrate1', help='substrate smile', type=str, required=True)
    parser.add_argument('--substrate2', help='substrate smile', type=str, required=True)
    parser.add_argument("--read_xyz", default=False, action="store_true" , help="To read or not to read?")
    args = parser.parse_args()
    main(args.substrate1,args.substrate2)
