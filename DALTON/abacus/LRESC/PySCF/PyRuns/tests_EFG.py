#import os 
#os.system("export PYTHONPATH=/home/juanjoaucar/Archivos_Juan/pyscf_new:$PYTHONPATH")
import pyscf
from pyscf import gto,scf
import numpy
#import pyscf.prop.efg.rhf as efg
#from pyscf.prop.efg import rhf as efg

#from properties.pyscf.prop.efg import rhf as efg

def _get_quadrupole_integrals(mol, atm_id):
    nao = mol.nao
    with mol.with_rinv_origin(mol.atom_coord(atm_id)):
        # Compute the integrals of quadrupole operator
        # (3 \vec{r} \vec{r} - r^2) / r^5
        ipipv = mol.intor('int1e_ipiprinv', 9).reshape(3,3,nao,nao)
        ipvip = mol.intor('int1e_iprinvip', 9).reshape(3,3,nao,nao)
        h1ao = ipipv + ipvip  # (nabla i | r/r^3 | j)
        h1ao = h1ao + h1ao.transpose(0,1,3,2)

#   Transformación de la traza que hace el pyscf (da lo mismo que Dalton)
#    coords = mol.atom_coord(atm_id).reshape(1, 3)
#    print('coordenadas',coords)
#    ao = mol.eval_gto('GTOval', coords)
#    print("ao",ao.shape,ao)
#    fc = 4*numpy.pi/3 * numpy.einsum('ip,iq->pq', ao, ao)
#    h1ao[0,0] += fc
#    h1ao[1,1] += fc
#    h1ao[2,2] += fc


    # Transformo la traza como hace Dalton
    qxx=numpy.copy(h1ao[0,0])
    qyy=numpy.copy(h1ao[1,1])
    qzz=numpy.copy(h1ao[2,2])

    h1ao[0,0] = (2*qxx-qyy-qzz)/3
    h1ao[1,1] = (2*qyy-qxx-qzz)/3
    h1ao[2,2] = (2*qzz-qxx-qyy)/3
    

    return h1ao

def _get_new_integrals(mol, atm_id, zeta = 9999999999):
    nao = mol.nao
#    with mol.with_rinv_origin(mol.atom_coord(atm_id)):
    with mol.with_rinv_zeta(zeta), mol.with_rinv_origin((mol.atom_coord(atm_id))):
        # Compute the integrals of quadrupole operator
        # (3 \vec{r} \vec{r} - r^2) / r^5
        c = mol.cart2sph_coeff()
        ipiprinvipip_cart = mol.intor('int1e_ipiprinvipip', 81).reshape(3,3,3,3,nao,nao)
        ipipiprinvip_cart = mol.intor('int1e_ipipiprinvip', 81).reshape(3,3,3,3,nao,nao)
        ipipipiprinv_cart = mol.intor('int1e_ipipipiprinv', 81).reshape(3,3,3,3,nao,nao)

        #Sphericals
#        nao_sph = mol_h2o.intor('int1e_ovlp_sph').shape[0]
        nao_sph = mol.intor('int1e_ovlp_sph').shape[0]
        ipiprinvipip = mol.intor('int1e_ipiprinvipip_sph', 81).reshape(3,3,3,3,nao_sph,nao_sph)
        ipipiprinvip = mol.intor('int1e_ipipiprinvip_sph', 81).reshape(3,3,3,3,nao_sph,nao_sph)
        ipipipiprinv = mol.intor('int1e_ipipipiprinv_sph', 81).reshape(3,3,3,3,nao_sph,nao_sph)

        #Chequeado para el elemento zzzz (2,2,2,2,:,:) con funciones s y con funciones p:
        core = -ipipiprinvip.transpose(3,1,2,0,5,4)-ipiprinvipip.transpose(0,3,1,2,4,5)-ipipiprinvip.transpose(1,3,2,0,5,4)-ipiprinvipip.transpose(3,0,1,2,4,5)-ipipiprinvip.transpose(1,2,3,0,5,4)-ipipipiprinv.transpose(1,2,3,0,5,4)-ipipiprinvip.transpose(0,1,2,3,4,5)-ipiprinvipip.transpose(0,1,3,2,4,5)

        #Nuevos índices (da lo mismo - chequeado)
#        core = -ipipiprinvip.transpose(3,1,0,2,5,4)-ipiprinvipip.transpose(1,2,0,3,4,5)-ipipiprinvip.transpose(1,3,0,2,5,4)-ipiprinvipip.transpose(2,1,0,3,4,5)-ipipiprinvip.transpose(1,0,3,2,5,4)-ipipipiprinv.transpose(2,1,0,3,5,4)-ipipiprinvip.transpose(2,1,0,3,4,5)-ipiprinvipip.transpose(1,0,2,3,4,5)

#        h1ao = ipipv + ipvip  # (nabla i | r/r^3 | j)
#        h1ao = h1ao + h1ao.transpose(0,1,3,2)

    return core
#    return h1ao

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
    factor=1.88972612499359 #Tomado de Dalton
    distance_bohr=factor*distance
    factorop=1/(3*distance_bohr*distance_bohr)
    coord_1=numpy.array([0,0,distance_bohr])
    coord_2=numpy.array([0,0,-distance_bohr])
    coord_3=numpy.array([0,distance_bohr,0])
    coord_4=numpy.array([0,-distance_bohr,0])
    coord_5=numpy.array([distance_bohr,0,0])
    coord_6=numpy.array([-distance_bohr,0,0])

    return numpy.array([coord_1,coord_2,coord_3,coord_4,coord_5,coord_6]),factorop



#mol = gto.M(atom = 'He 0 0 0', basis = 'ccpvdz')
#mol = gto.M(atom = 'He 0 0 0; He 0 0 2', basis = 'ccpvdz')

#with open('molecula.xyz', 'r') as f:
#    base = f.readlines()[-1]
base = "dyall_cv2z"

#mol_h2o = gto.M(atom="molecula.xyz",basis = base[:-1])
mol_h2o = gto.M(atom="molecula.xyz",basis = base)
mol_h2o, ctr_coeff = mol_h2o.to_uncontracted_cartesian_basis()

#Tests transformación a esféricas
#print("atomic orbitals:",mol_h2o.nao, mol_h2o.nao_cart(),mol_h2o.nao_nr(cart=None))
#c = mol_h2o.cart2sph_coeff()
#print("test_sph1:",mol_h2o.intor('int1e_ovlp_sph').shape)
#print("test_sph2",c.T.dot(mol_h2o.intor('int1e_ovlp_cart')).dot(c).shape)
#print("test_cart:",mol_h2o.intor('int1e_ovlp_cart').shape)
#naos_sph = mol_h2o.intor('int1e_ovlp_sph').shape[0]
#print("naos_sph",naos_sph)
#print("nbas",mol_h2o.nbas)


#rhf_h2o = scf.RHF(mol_h2o)
#e_h2o = rhf_h2o.kernel()


#pyscf.prop.efg.rhf.EFG(method, efg_nuc=None)¶
#pyscf.prop.efg.rhf.kernel(method, efg_nuc=None)¶

#atm_id=0
#integrals=efg._get_quadrupole_integrals(mol_h2o, atm_id)
#integrals=_get_quadrupole_integrals(mol_h2o, atm_id)
#print(integrals.shape)
#print("zz integrals:")
#print(integrals[2,2,:,:])


#integrals_new=_get_new_integrals(mol_h2o, atm_id)
#print(integrals_new.shape)
#DERIVADA Z Z Z Z
#print(integrals_new[2,2,2,2,:,:])
#print(integrals_new[2,2,2,2,:,:].shape)

orbitales = numpy.array([mol_h2o.nao])
naos_sph = mol_h2o.intor('int1e_ovlp_sph').shape[0]
nro_operadores=numpy.array([3*mol_h2o.natm])
print("naos_cart:",orbitales)
print("naos_sph:",naos_sph)
print("nro de operadores a grabar: ",nro_operadores)
with open("shi.bin", "wb") as file:
    numpy.array(nro_operadores[0], dtype=numpy.int32).tofile(file)
#    numpy.array(orbitales[0], dtype=numpy.int32).tofile(file)
    numpy.array(naos_sph, dtype=numpy.int32).tofile(file)

for i in range(mol_h2o.natm):
    titulos=['newopx'+"{:02.0f}".format(i+1), 'newopy'+"{:02.0f}".format(i+1), 'newopz'+"{:02.0f}".format(i+1)]
    print(titulos)
    valorzeta=EXPONENTEMANUAL #Para variar en script
#    valorzeta=10E+11 #Modelo puntual
#    valorzeta=1E+05
#    valorzeta=5.3549302238E+08 # for F
#    valorzeta=2.4319909606E+08 #for Br
#    valorzeta=1.8452373543E+08 #for I
#   Agrego el signo menos debido a las p
    nuevo_x=-(_get_GiZZj(mol_h2o,i,1,2, zeta=valorzeta)-_get_GiZZj(mol_h2o,i,2,1, zeta=valorzeta))
    nuevo_y=-(_get_GiZZj(mol_h2o,i,2,0, zeta=valorzeta)-_get_GiZZj(mol_h2o,i,0,2, zeta=valorzeta))
    nuevo_z=-(_get_GiZZj(mol_h2o,i,0,1, zeta=valorzeta)-_get_GiZZj(mol_h2o,i,1,0, zeta=valorzeta))

#    with open("test_c.bin", "wb") as file:
#        numpy.array(nro_operadores[0], dtype=numpy.int32).tofile(file)
#        numpy.array(orbitales[0], dtype=numpy.int32).tofile(file)

#    for j in range(nro_operadores[0]):
#    somestring1 = titulos[i]
    with open("shi.bin", "ab") as file:
        file.write(titulos[0].encode('ascii'))
        numpy.array(nuevo_x, dtype=numpy.float64).transpose(1,0).tofile(file)
        print("grabando en shi.bin el operador ", titulos[0])
        file.write(titulos[1].encode('ascii'))
        numpy.array(nuevo_y, dtype=numpy.float64).transpose(1,0).tofile(file)
        file.write(titulos[2].encode('ascii'))
        numpy.array(nuevo_z, dtype=numpy.float64).transpose(1,0).tofile(file)

"""
#SECCION PCNQM
distance=0.0001 #In Angstroms
charges_coord,factorPCNQM=_PCNQM_charges(distance)

PCNQM_x=[0,0,0,0,0,0]
PCNQM_y=[0,0,0,0,0,0]
PCNQM_z=[0,0,0,0,0,0]
for i in range(6):
    PCNQM_x[i],PCNQM_y[i],PCNQM_z[i]=_get_newop_SO_S_qzz_PCNQM(mol_h2o,charges_coord,i, zeta = valorzeta)
#QUITO EL FACTOR (se lo agrega luego de la LR):
factorPCNQM=1.0

newop_PCNQM_x=factorPCNQM*(-2*PCNQM_x[0]-2*PCNQM_x[1]+PCNQM_x[2]+PCNQM_x[3]+PCNQM_x[4]+PCNQM_x[5])
newop_PCNQM_y=factorPCNQM*(-2*PCNQM_y[0]-2*PCNQM_y[1]+PCNQM_y[2]+PCNQM_y[3]+PCNQM_y[4]+PCNQM_y[5])
newop_PCNQM_z=factorPCNQM*(-2*PCNQM_z[0]-2*PCNQM_z[1]+PCNQM_z[2]+PCNQM_z[3]+PCNQM_z[4]+PCNQM_z[5])


with open("test_c.bin", "ab") as file:
    file.write('PCNQM  x'.encode('ascii'))
    numpy.array(newop_PCNQM_x, dtype=numpy.float64).transpose(1,0).tofile(file)
    file.write('PCNQM  y'.encode('ascii'))
    numpy.array(newop_PCNQM_y, dtype=numpy.float64).transpose(1,0).tofile(file)
    file.write('PCNQM  z'.encode('ascii'))
    numpy.array(newop_PCNQM_z, dtype=numpy.float64).transpose(1,0).tofile(file)
"""
