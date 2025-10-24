# ######################################################
#
# Writen by Juan J. Aucar. August 2025.
#
# Generalization of the FC mechanism. See Eq. 25 at
# https://doi.org/10.1063/5.0264596
#
# Feel free to contact me at juanaucar@gmail.com
#
# ######################################################
import pyscf
from pyscf import gto,scf
import numpy

def _get_DIFC_I(mol, atm_id, cartesian = False ):
    """Calculate (r_i . nabla_j) / (r^3) integrals.

    :param mol: The molecular object.
    :type mol: pyscf.gto.Mole
    :param atm_id: The atom index.
    :type atm_id: int
    :param cartesian: Whether to use Cartesian gaussians, defaults to False
    :type cartesian: bool, optional
    """
    nao_sph = mol.intor('int1e_ovlp_sph').shape[0]
    nao_cart = mol.intor('int1e_ovlp_cart').shape[0]

    with mol.with_rinv_origin((mol.atom_coord(atm_id))):
        if (cartesian):
#	    Cartesians
            iprinvip = mol.intor('int1e_iprinvip_cart', 9).reshape(3,3,nao_cart,nao_cart)
            ipiprinv = mol.intor('int1e_ipiprinv_cart', 9).reshape(3,3,nao_cart,nao_cart)
            rinvipip=numpy.transpose(ipiprinv,axes=[0,1,3,2]).conjugate()
        else :
#	    Sphericals
            iprinvip = mol.intor('int1e_iprinvip_sph', 9).reshape(3,3,nao_sph,nao_sph)
            ipiprinv = mol.intor('int1e_ipiprinv_sph', 9).reshape(3,3,nao_sph,nao_sph)
            rinvipip=numpy.transpose(ipiprinv,axes=[0,1,3,2]).conjugate()

    integrals=iprinvip+rinvipip #ipiprinv
    return integrals

def _get_DIFC_0(mol, atm_id, cartesian = False ):
    """Calculate (r . nabla) / (r^3) integrals.

    :param mol: The molecular object.
    :type mol: pyscf.gto.Mole
    :param atm_id: The atom index.
    :type atm_id: int
    :param cartesian: Whether to use Cartesian gaussians, defaults to False
    :type cartesian: bool, optional
    :return: The (r . nabla) / (r^3) integrals.
    :rtype: numpy.ndarray
    """
    DIFC_I = _get_DIFC_I(mol, atm_id, cartesian)
    return DIFC_I[0,0,:,:]+DIFC_I[1,1,:,:]+DIFC_I[2,2,:,:]

def _get_DIFC_II(mol, atm_id, cartesian = False ):
    """Calculate (r_i r_j)/(r^5) r . nabla integrals.

    :param mol: The molecular object.
    :type mol: pyscf.gto.Mole
    :param atm_id: The atom index.
    :type atm_id: int
    :param cartesian: Whether to use Cartesian gaussians, defaults to False
    :type cartesian: bool, optional
    :return: The (r_i r_j)/(r^5) r . nabla integrals.
    :rtype: numpy.ndarray
    """
    nao_sph = mol.intor('int1e_ovlp_sph').shape[0]
    nao_cart = mol.intor('int1e_ovlp_cart').shape[0]

    with mol.with_rinv_origin((mol.atom_coord(atm_id))):
        if (cartesian):
#	    Cartesians
            iprinviprip = mol_system.intor('int1e_iprinviprip_cart', 81).reshape(3,3,3,3,nao_cart,nao_cart)
            rinvipiprip = mol_system.intor('int1e_rinvipiprip_cart', 81).reshape(3,3,3,3,nao_cart,nao_cart)
            ipiprinvrip = mol_system.intor('int1e_ipiprinvrip_cart', 81).reshape(3,3,3,3,nao_cart,nao_cart)
            iprinvip = mol.intor('int1e_iprinvip_cart', 9).reshape(3,3,nao_cart,nao_cart)
            ipiprinv = mol.intor('int1e_ipiprinv_cart', 9).reshape(3,3,nao_cart,nao_cart)
            rinvipip=numpy.transpose(ipiprinv,axes=[0,1,3,2]).conjugate()
            integrals=numpy.zeros([3,3,nao_cart,nao_cart])
        else :
#	    Sphericals
            iprinviprip = mol_system.intor('int1e_iprinviprip_sph', 81).reshape(3,3,3,3,nao_sph,nao_sph)
            rinvipiprip = mol_system.intor('int1e_rinvipiprip_sph', 81).reshape(3,3,3,3,nao_sph,nao_sph)
            ipiprinvrip = mol_system.intor('int1e_ipiprinvrip_sph', 81).reshape(3,3,3,3,nao_sph,nao_sph)
            iprinvip = mol.intor('int1e_iprinvip_sph', 9).reshape(3,3,nao_sph,nao_sph)
            ipiprinv = mol.intor('int1e_ipiprinv_sph', 9).reshape(3,3,nao_sph,nao_sph)
            rinvipip=numpy.transpose(ipiprinv,axes=[0,1,3,2]).conjugate()
            integrals=numpy.zeros([3,3,nao_sph,nao_sph])

    for k in range(3):
        integrals = integrals + iprinviprip[:,:,k,k,:,:] + numpy.transpose(iprinviprip[:,:,k,k,:,:],axes=[1,0,2,3]) + numpy.transpose(rinvipiprip[:,:,k,k,:,:],axes=[1,0,2,3]) + numpy.transpose(ipiprinvrip[:,:,k,k,:,:],axes=[1,0,2,3]) #+ iprinvip[k,k,:,:] + rinvipip[k,k,:,:]

    for k in range(3):
        integrals[0,0,:,:] = integrals[0,0,:,:] + iprinvip[k,k,:,:] + rinvipip[k,k,:,:]
        integrals[1,1,:,:] = integrals[1,1,:,:] + iprinvip[k,k,:,:] + rinvipip[k,k,:,:]
        integrals[2,2,:,:] = integrals[2,2,:,:] + iprinvip[k,k,:,:] + rinvipip[k,k,:,:]

    return integrals/3


# ######################################################
# Main code starts here
# ######################################################


base = "dyall_cv4z"

mol_system = gto.M(atom="molecula.xyz",basis = base)
mol_system, ctr_coeff = mol_system.to_uncontracted_cartesian_basis()

naos_sph = mol_system.intor('int1e_ovlp_sph').shape[0]
nao_cart = mol_system.intor('int1e_ovlp_cart').shape[0]
operators_N=numpy.array([19*mol_system.natm])
print("naos_cart:",nao_cart)
print("naos_sph:",naos_sph)
print("Number of operators: ",operators_N)

cartesian = False
if (cartesian):
   print("Cartesian basis set used at PySCF")
else:
   print("Spherical basis set used at PySCF")


with open("shi.bin", "wb") as file:
    numpy.array(operators_N[0], dtype=numpy.int32).tofile(file)
    if (cartesian):
       numpy.array(nao_cart, dtype=numpy.int32).tofile(file)
    else:
       numpy.array(naos_sph, dtype=numpy.int32).tofile(file)


# Iterate over atoms to obtain the integrals
for i in range(mol_system.natm):
    DIFC_0=_get_DIFC_0(mol_system, atm_id = i, cartesian = cartesian)
    titulos1=['FC1  '+"{:03.0f}".format(i+1)]

    DIFC_I=_get_DIFC_I(mol_system, atm_id = i, cartesian = cartesian)
    titulos2_row1=['FC2xx'+"{:03.0f}".format(i+1), 'FC2xy'+"{:03.0f}".format(i+1), 'FC2xz'+"{:03.0f}".format(i+1)]
    titulos2_row2=['FC2yx'+"{:03.0f}".format(i+1), 'FC2yy'+"{:03.0f}".format(i+1), 'FC2yz'+"{:03.0f}".format(i+1)]
    titulos2_row3=['FC2zx'+"{:03.0f}".format(i+1), 'FC2zy'+"{:03.0f}".format(i+1), 'FC2zz'+"{:03.0f}".format(i+1)]

    DIFC_II=_get_DIFC_II(mol_system, atm_id = i, cartesian = cartesian)
    titulos3_row1=['FC3xx'+"{:03.0f}".format(i+1), 'FC3xy'+"{:03.0f}".format(i+1), 'FC3xz'+"{:03.0f}".format(i+1)]
    titulos3_row2=['FC3yx'+"{:03.0f}".format(i+1), 'FC3yy'+"{:03.0f}".format(i+1), 'FC3yz'+"{:03.0f}".format(i+1)]
    titulos3_row3=['FC3zx'+"{:03.0f}".format(i+1), 'FC3zy'+"{:03.0f}".format(i+1), 'FC3zz'+"{:03.0f}".format(i+1)]

    with open("shi.bin", "ab") as file:
        file.write(titulos1[0].encode('ascii'))
        numpy.array(DIFC_0[:,:], dtype=numpy.float64).transpose(1,0).tofile(file)

        file.write(titulos2_row1[0].encode('ascii'))
        numpy.array(DIFC_I[0,0,:,:], dtype=numpy.float64).transpose(1,0).tofile(file)
        file.write(titulos2_row1[1].encode('ascii'))
        numpy.array(DIFC_I[0,1,:,:], dtype=numpy.float64).transpose(1,0).tofile(file)
        file.write(titulos2_row1[2].encode('ascii'))
        numpy.array(DIFC_I[0,2,:,:], dtype=numpy.float64).transpose(1,0).tofile(file)
        file.write(titulos2_row2[0].encode('ascii'))
        numpy.array(DIFC_I[1,0,:,:], dtype=numpy.float64).transpose(1,0).tofile(file)
        file.write(titulos2_row2[1].encode('ascii'))
        numpy.array(DIFC_I[1,1,:,:], dtype=numpy.float64).transpose(1,0).tofile(file)
        file.write(titulos2_row2[2].encode('ascii'))
        numpy.array(DIFC_I[1,2,:,:], dtype=numpy.float64).transpose(1,0).tofile(file)
        file.write(titulos2_row3[0].encode('ascii'))
        numpy.array(DIFC_I[2,0,:,:], dtype=numpy.float64).transpose(1,0).tofile(file)
        file.write(titulos2_row3[1].encode('ascii'))
        numpy.array(DIFC_I[2,1,:,:], dtype=numpy.float64).transpose(1,0).tofile(file)
        file.write(titulos2_row3[2].encode('ascii'))
        numpy.array(DIFC_I[2,2,:,:], dtype=numpy.float64).transpose(1,0).tofile(file)

        file.write(titulos3_row1[0].encode('ascii'))
        numpy.array(DIFC_II[0,0,:,:], dtype=numpy.float64).transpose(1,0).tofile(file)
        file.write(titulos3_row1[1].encode('ascii'))
        numpy.array(DIFC_II[0,1,:,:], dtype=numpy.float64).transpose(1,0).tofile(file)
        file.write(titulos3_row1[2].encode('ascii'))
        numpy.array(DIFC_II[0,2,:,:], dtype=numpy.float64).transpose(1,0).tofile(file)
        file.write(titulos3_row2[0].encode('ascii'))
        numpy.array(DIFC_II[1,0,:,:], dtype=numpy.float64).transpose(1,0).tofile(file)
        file.write(titulos3_row2[1].encode('ascii'))
        numpy.array(DIFC_II[1,1,:,:], dtype=numpy.float64).transpose(1,0).tofile(file)
        file.write(titulos3_row2[2].encode('ascii'))
        numpy.array(DIFC_II[1,2,:,:], dtype=numpy.float64).transpose(1,0).tofile(file)
        file.write(titulos3_row3[0].encode('ascii'))
        numpy.array(DIFC_II[2,0,:,:], dtype=numpy.float64).transpose(1,0).tofile(file)
        file.write(titulos3_row3[1].encode('ascii'))
        numpy.array(DIFC_II[2,1,:,:], dtype=numpy.float64).transpose(1,0).tofile(file)
        file.write(titulos3_row3[2].encode('ascii'))
        numpy.array(DIFC_II[2,2,:,:], dtype=numpy.float64).transpose(1,0).tofile(file)
