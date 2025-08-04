# ######################################################
#
# Writen by Juan J. Aucar. August 2025.
#
# Results for the SOS correction were first published at
# https://doi.org/10.1063/5.0124701
#
# This file contains some tests that may be of interest
# to the user. 
#
# Feel free to contact me at juanaucar@gmail.com
#
# ######################################################
import pyscf
from pyscf import gto,scf
import numpy

def _get_quadrupole_integrals(mol, atm_id):
    nao = mol.nao
    with mol.with_rinv_origin(mol.atom_coord(atm_id)):
        # Compute the integrals of quadrupole operator
        # (3 \vec{r} \vec{r} - r^2) / r^5
        ipipv = mol.intor('int1e_ipiprinv', 9).reshape(3,3,nao,nao)
        ipvip = mol.intor('int1e_iprinvip', 9).reshape(3,3,nao,nao)
        h1ao = ipipv + ipvip  # (nabla i | r/r^3 | j)
        h1ao = h1ao + h1ao.transpose(0,1,3,2)

#    Transformation to null trace tensor as PySCF does
#    coords = mol.atom_coord(atm_id).reshape(1, 3)
#    print('coordenadas',coords)
#    ao = mol.eval_gto('GTOval', coords)
#    print("ao",ao.shape,ao)
#    fc = 4*numpy.pi/3 * numpy.einsum('ip,iq->pq', ao, ao)
#    h1ao[0,0] += fc
#    h1ao[1,1] += fc
#    h1ao[2,2] += fc

    #    Transformation to null trace tensor as DALTON does
    qxx=numpy.copy(h1ao[0,0])
    qyy=numpy.copy(h1ao[1,1])
    qzz=numpy.copy(h1ao[2,2])

    h1ao[0,0] = (2*qxx-qyy-qzz)/3
    h1ao[1,1] = (2*qyy-qxx-qzz)/3
    h1ao[2,2] = (2*qzz-qxx-qyy)/3
    
    return h1ao

# See function defined with the same name at EFG_SOS_sph.py file.
def _get_new_integrals(mol, atm_id, zeta = 9999999999):
    nao = mol.nao
    with mol.with_rinv_zeta(zeta), mol.with_rinv_origin((mol.atom_coord(atm_id))):
        # Compute the integrals of quadrupole operator
        # (3 \vec{r} \vec{r} - r^2) / r^5
        c = mol.cart2sph_coeff()
        ipiprinvipip_cart = mol.intor('int1e_ipiprinvipip', 81).reshape(3,3,3,3,nao,nao)
        ipipiprinvip_cart = mol.intor('int1e_ipipiprinvip', 81).reshape(3,3,3,3,nao,nao)
        ipipipiprinv_cart = mol.intor('int1e_ipipipiprinv', 81).reshape(3,3,3,3,nao,nao)

        #Sphericals
        nao_sph = mol.intor('int1e_ovlp_sph').shape[0]
        ipiprinvipip = mol.intor('int1e_ipiprinvipip_sph', 81).reshape(3,3,3,3,nao_sph,nao_sph)
        ipipiprinvip = mol.intor('int1e_ipipiprinvip_sph', 81).reshape(3,3,3,3,nao_sph,nao_sph)
        ipipipiprinv = mol.intor('int1e_ipipipiprinv_sph', 81).reshape(3,3,3,3,nao_sph,nao_sph)

        integrals = -ipipiprinvip.transpose(3,1,2,0,5,4)-ipiprinvipip.transpose(0,3,1,2,4,5)-ipipiprinvip.transpose(1,3,2,0,5,4)-ipiprinvipip.transpose(3,0,1,2,4,5)-ipipiprinvip.transpose(1,2,3,0,5,4)-ipipipiprinv.transpose(1,2,3,0,5,4)-ipipiprinvip.transpose(0,1,2,3,4,5)-ipiprinvipip.transpose(0,1,3,2,4,5)


#        h1ao = ipipv + ipvip  # (nabla i | r/r^3 | j)
#        h1ao = h1ao + h1ao.transpose(0,1,3,2)

    return integrals

# To do some TESTS related to the PCNQM method
def _get_new_integrals_PCNQM(mol, charges_coord,atm_id, zeta = 9999999999):
    nao = mol.nao
    with mol.with_rinv_zeta(zeta), mol.with_rinv_origin(charges_coord[atm_id]):
        ipiprinv = mol.intor('int1e_ipiprinv', 9).reshape(3,3,nao,nao)
        iprinvip = mol.intor('int1e_iprinvip', 9).reshape(3,3,nao,nao)
        core = -iprinvip.transpose(0,1,2,3)-ipiprinv.transpose(0,1,3,2)

    return core

def _get_GiZZj(mol,atm_id,i,j, zeta = 9999999999):
    #Gets the transformed matrix
    matriz_xx=_get_new_integrals(mol,atm_id,zeta)[i,0,0,j,:,:]
    matriz_yy=_get_new_integrals(mol,atm_id,zeta)[i,1,1,j,:,:]
    matriz_zz=_get_new_integrals(mol,atm_id,zeta)[i,2,2,j,:,:]
    return 1/3*(2*matriz_zz-matriz_xx-matriz_yy)


def _get_newop_SO_S_qzz_PCNQM(mol,charges_coord,atm_id, zeta = 9999999999):
    #Gets the untransformed matrix for qzz SO_S  (no need to transform)
    matrix_x=-(_get_new_integrals_PCNQM(mol,charges_coord,atm_id,zeta)[1,2,:,:]
             -_get_new_integrals_PCNQM(mol,charges_coord,atm_id,zeta)[2,1,:,:])
    matrix_y=-(_get_new_integrals_PCNQM(mol,charges_coord,atm_id,zeta)[2,0,:,:]
             -_get_new_integrals_PCNQM(mol,charges_coord,atm_id,zeta)[0,2,:,:])
    matrix_z=-(_get_new_integrals_PCNQM(mol,charges_coord,atm_id,zeta)[0,1,:,:]
             -_get_new_integrals_PCNQM(mol,charges_coord,atm_id,zeta)[1,0,:,:])
    return matrix_x,matrix_y,matrix_z


def _PCNQM_charges(distance):
    #Input distance: Angstrom
    #Output distance: a.u. (bohr)
    factor=1.88972612499359 # Taken from DALTON
    distance_bohr=factor*distance
    factorop=1/(3*distance_bohr*distance_bohr)
    coord_1=numpy.array([0,0,distance_bohr])
    coord_2=numpy.array([0,0,-distance_bohr])
    coord_3=numpy.array([0,distance_bohr,0])
    coord_4=numpy.array([0,-distance_bohr,0])
    coord_5=numpy.array([distance_bohr,0,0])
    coord_6=numpy.array([-distance_bohr,0,0])

    return numpy.array([coord_1,coord_2,coord_3,coord_4,coord_5,coord_6]),factorop


# ######################################################
# Test code starts here
# ######################################################


#mol = gto.M(atom = 'He 0 0 0', basis = 'ccpvdz') # This can be entered in the command line
base = "dyall_cv2z"

mol_h2o = gto.M(atom="molecula.xyz",basis = base)
mol_h2o, ctr_coeff = mol_h2o.to_uncontracted_cartesian_basis()

#rhf_h2o = scf.RHF(mol_h2o)
#e_h2o = rhf_h2o.kernel()

#pyscf.prop.efg.rhf.EFG(method, efg_nuc=None)
#pyscf.prop.efg.rhf.kernel(method, efg_nuc=None)

#atm_id=0
#integrals=efg._get_quadrupole_integrals(mol_h2o, atm_id)
#integrals=_get_quadrupole_integrals(mol_h2o, atm_id)
#print(integrals.shape)
#print("zz integrals:")
#print(integrals[2,2,:,:])

#integrals_new=_get_new_integrals(mol_h2o, atm_id)
#print(integrals_new.shape)
#Derivative with respect to Z Z Z Z
#print(integrals_new[2,2,2,2,:,:])
#print(integrals_new[2,2,2,2,:,:].shape)

orbitals = numpy.array([mol_h2o.nao])
naos_sph = mol_h2o.intor('int1e_ovlp_sph').shape[0]
nro_operadores=numpy.array([3*mol_h2o.natm])
print("naos_cart:",orbitals)
print("naos_sph:",naos_sph)
print("Number of operators: ",nro_operadores)


"""
# PCNQM method
distance=0.0001 #In Angstroms
charges_coord,factorPCNQM=_PCNQM_charges(distance)

PCNQM_x=[0,0,0,0,0,0]
PCNQM_y=[0,0,0,0,0,0]
PCNQM_z=[0,0,0,0,0,0]
for i in range(6):
    PCNQM_x[i],PCNQM_y[i],PCNQM_z[i]=_get_newop_SO_S_qzz_PCNQM(mol_h2o,charges_coord,i, zeta = valorzeta)
factorPCNQM=1.0

newop_PCNQM_x=factorPCNQM*(-2*PCNQM_x[0]-2*PCNQM_x[1]+PCNQM_x[2]+PCNQM_x[3]+PCNQM_x[4]+PCNQM_x[5])
newop_PCNQM_y=factorPCNQM*(-2*PCNQM_y[0]-2*PCNQM_y[1]+PCNQM_y[2]+PCNQM_y[3]+PCNQM_y[4]+PCNQM_y[5])
newop_PCNQM_z=factorPCNQM*(-2*PCNQM_z[0]-2*PCNQM_z[1]+PCNQM_z[2]+PCNQM_z[3]+PCNQM_z[4]+PCNQM_z[5])

"""
