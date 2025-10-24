# ######################################################
#
# Writen by Juan J. Aucar. August 2025.
#
# Results for the SOS correction were first published at
# https://doi.org/10.1063/5.0124701
#
# Feel free to contact me at juanaucar@gmail.com
#
# ######################################################
from pyscf import gto,scf
import numpy

def _get_new_integrals(mol, atm_id, zeta = 9999999999, cartesian = False):
    """Get integrals for a given atom in the molecule.

    :param mol: The molecule object.
    :type mol: pyscf.gto.Mole
    :param atm_id: The atom ID for which to compute the integrals.
    :type atm_id: int
    :param zeta: The zeta parameter for the integral calculation, defaults to 9999999999
    :type zeta: int, optional
    :param cartesian: Whether to use Cartesian gaussians, defaults to False
    :type cartesian: bool, optional
    """
    #Get the dimension of the atomic basis
    nao_cart = mol.intor('int1e_ovlp_cart').shape[0]
    nao_sph = mol.intor('int1e_ovlp_sph').shape[0]
    with mol.with_rinv_zeta(zeta), mol.with_rinv_origin((mol.atom_coord(atm_id))):
        if (cartesian):
 	    #Cartesians
            ipiprinvipip = mol.intor('int1e_ipiprinvipip', 81).reshape(3,3,3,3,nao_cart,nao_cart)
            ipipiprinvip = mol.intor('int1e_ipipiprinvip', 81).reshape(3,3,3,3,nao_cart,nao_cart)
            ipipipiprinv = mol.intor('int1e_ipipipiprinv', 81).reshape(3,3,3,3,nao_cart,nao_cart)
        else:
	    #Sphericals
            ipiprinvipip = mol.intor('int1e_ipiprinvipip_sph', 81).reshape(3,3,3,3,nao_sph,nao_sph)
            ipipiprinvip = mol.intor('int1e_ipipiprinvip_sph', 81).reshape(3,3,3,3,nao_sph,nao_sph)
            ipipipiprinv = mol.intor('int1e_ipipipiprinv_sph', 81).reshape(3,3,3,3,nao_sph,nao_sph)

        integrals = -ipipiprinvip.transpose(3,1,2,0,5,4)-ipiprinvipip.transpose(0,3,1,2,4,5)-ipipiprinvip.transpose(1,3,2,0,5,4)-ipiprinvipip.transpose(3,0,1,2,4,5)-ipipiprinvip.transpose(1,2,3,0,5,4)-ipipipiprinv.transpose(1,2,3,0,5,4)-ipipiprinvip.transpose(0,1,2,3,4,5)-ipiprinvipip.transpose(0,1,3,2,4,5)

        #Transpose (gives same result)
#        integrals = -ipipiprinvip.transpose(3,1,0,2,5,4)-ipiprinvipip.transpose(1,2,0,3,4,5)-ipipiprinvip.transpose(1,3,0,2,5,4)-ipiprinvipip.transpose(2,1,0,3,4,5)-ipipiprinvip.transpose(1,0,3,2,5,4)-ipipipiprinv.transpose(2,1,0,3,5,4)-ipipiprinvip.transpose(2,1,0,3,4,5)-ipiprinvipip.transpose(1,0,2,3,4,5)


    return integrals

def _get_GiZZj(mol,atm_id,i,j, zeta = 9999999999, cartesian = False):
    """Generate the integrals for a given atom in the molecule.

    :param mol: The molecule object.
    :type mol: pyscf.gto.Mole
    :param atm_id: The atom ID for which to compute the integrals.
    :type atm_id: int
    :param i: The index i for the integral.
    :type i: int
    :param j: The index j for the integral.
    :type j: int
    :param zeta: The zeta parameter for the integral calculation, defaults to 9999999999
    :type zeta: int, optional
    :param cartesian: Whether to use Cartesian gaussians, defaults to False
    :type cartesian: bool, optional
    :return: The G^i_{zz}^j integral.
    :rtype: numpy.ndarray
    """
    #Gets the transformed matrix
    matriz_xx=_get_new_integrals(mol,atm_id,zeta,cartesian)[i,0,0,j,:,:]
    matriz_yy=_get_new_integrals(mol,atm_id,zeta,cartesian)[i,1,1,j,:,:]
    matriz_zz=_get_new_integrals(mol,atm_id,zeta,cartesian)[i,2,2,j,:,:]
    return 1/3*(2*matriz_zz-matriz_xx-matriz_yy)


# ######################################################
# Main code starts here
# ######################################################

base = "dyall_cv4z"

mol_system = gto.M(atom="molecula.xyz",basis = base)
mol_system, ctr_coeff = mol_system.to_uncontracted_cartesian_basis()

naos_sph = mol_system.intor('int1e_ovlp_sph').shape[0]
nao_cart = mol_system.intor('int1e_ovlp_cart').shape[0]
operators_N=numpy.array([3*mol_system.natm])
print("naos_cart:",nao_cart)
print("naos_sph:",naos_sph)
print("Number of operators to save: ",operators_N)
cartesian = False
if (cartesian):
   print("Cartesian basis set used at PySCF")
else:
   print("Spherical basis set used at PySCF")

with open("efg.bin", "wb") as file:
    numpy.array(operators_N[0], dtype=numpy.int32).tofile(file)
    if (cartesian):
       numpy.array(nao_cart, dtype=numpy.int32).tofile(file)
    else:
       numpy.array(naos_sph, dtype=numpy.int32).tofile(file)


for i in range(mol_system.natm):
    titulos=['newopx'+"{:02.0f}".format(i+1), 'newopy'+"{:02.0f}".format(i+1), 'newopz'+"{:02.0f}".format(i+1)]
    print(titulos)
    valorzeta=EXPONENTEMANUAL #Change this value to the proper zeta
#    valorzeta=10E+11 #Point Model


#   Minus sign added due to p
    new_x=-(_get_GiZZj(mol_system,i,1,2, zeta=valorzeta, cartesian = cartesian)-_get_GiZZj(mol_system,i,2,1, zeta=valorzeta, cartesian = cartesian))
    new_y=-(_get_GiZZj(mol_system,i,2,0, zeta=valorzeta, cartesian = cartesian)-_get_GiZZj(mol_system,i,0,2, zeta=valorzeta, cartesian = cartesian))
    new_z=-(_get_GiZZj(mol_system,i,0,1, zeta=valorzeta, cartesian = cartesian)-_get_GiZZj(mol_system,i,1,0, zeta=valorzeta, cartesian = cartesian))

    with open("efg.bin", "ab") as file:
        file.write(titulos[0].encode('ascii'))
        numpy.array(new_x, dtype=numpy.float64).transpose(1,0).tofile(file)
        print("Saving the operator in file efg.bin ", titulos[0])
        file.write(titulos[1].encode('ascii'))
        numpy.array(new_y, dtype=numpy.float64).transpose(1,0).tofile(file)
        file.write(titulos[2].encode('ascii'))
        numpy.array(new_z, dtype=numpy.float64).transpose(1,0).tofile(file)