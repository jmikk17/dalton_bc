import pyscf
from pyscf import gto,scf
import numpy

def _get_DIFC_I(mol, atm_id, cartesian = False ):
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

    integrals=iprinvip+ipiprinv
    return integrals


def _get_DIFC_II(mol, atm_id, cartesian = False ):
    nao_sph = mol.intor('int1e_ovlp_sph').shape[0]
    nao_cart = mol.intor('int1e_ovlp_cart').shape[0]

    with mol.with_rinv_origin((mol.atom_coord(atm_id))):
        if (cartesian):
#	    Cartesians
            iprinviprip = mol_h2o.intor('int1e_iprinviprip_cart', 81).reshape(3,3,3,3,nao_cart,nao_cart)
            rinvipiprip = mol_h2o.intor('int1e_rinvipiprip_cart', 81).reshape(3,3,3,3,nao_cart,nao_cart)
            ipiprinvrip = mol_h2o.intor('int1e_ipiprinvrip_cart', 81).reshape(3,3,3,3,nao_cart,nao_cart)
            iprinvip = mol.intor('int1e_iprinvip_cart', 9).reshape(3,3,nao_cart,nao_cart)
            ipiprinv = mol.intor('int1e_ipiprinv_cart', 9).reshape(3,3,nao_cart,nao_cart)
            rinvipip=numpy.transpose(ipiprinv,axes=[0,1,3,2]).conjugate()
            integrals=numpy.zeros([3,3,nao_cart,nao_cart])
        else :
#	    Sphericals
            iprinviprip = mol_h2o.intor('int1e_iprinviprip_sph', 81).reshape(3,3,3,3,nao_sph,nao_sph)
            rinvipiprip = mol_h2o.intor('int1e_rinvipiprip_sph', 81).reshape(3,3,3,3,nao_sph,nao_sph)
            ipiprinvrip = mol_h2o.intor('int1e_ipiprinvrip_sph', 81).reshape(3,3,3,3,nao_sph,nao_sph)
            iprinvip = mol.intor('int1e_iprinvip_sph', 9).reshape(3,3,nao_sph,nao_sph)
            ipiprinv = mol.intor('int1e_ipiprinv_sph', 9).reshape(3,3,nao_sph,nao_sph)
            rinvipip=numpy.transpose(ipiprinv,axes=[0,1,3,2]).conjugate()
            integrals=numpy.zeros([3,3,nao_sph,nao_sph])

    for k in range(3):
        integrals = integrals + iprinviprip[:,:,k,k,:,:] + numpy.transpose(iprinviprip[:,:,k,k,:,:],axes=[1,0,2,3]) + numpy.transpose(rinvipiprip[:,:,k,k,:,:],axes=[1,0,2,3]) + numpy.transpose(ipiprinvrip[:,:,k,k,:,:],axes=[1,0,2,3]) +iprinvip[k,k,:,:] + rinvipip[k,k,:,:]


    return integrals/3




base = "dyall_cv2z"

mol_h2o = gto.M(atom="molecula.xyz",basis = base)
mol_h2o, ctr_coeff = mol_h2o.to_uncontracted_cartesian_basis()

naos_sph = mol_h2o.intor('int1e_ovlp_sph').shape[0]
nao_cart = mol_h2o.intor('int1e_ovlp_cart').shape[0]
nro_operadores=numpy.array([18*mol_h2o.natm])
print("naos_cart:",nao_cart)
print("naos_sph:",naos_sph)
print("nro de operadores a grabar: ",nro_operadores)
cartesian = True
if (cartesian):
   print("Cartesian basis set used at PySCF")
else:
   print("Spherical basis set used at PySCF")


with open("shi.bin", "wb") as file:
    numpy.array(nro_operadores[0], dtype=numpy.int32).tofile(file)
    if (cartesian):
       numpy.array(nao_cart, dtype=numpy.int32).tofile(file)
    else:
       numpy.array(naos_sph, dtype=numpy.int32).tofile(file)



for i in range(mol_h2o.natm):
    DIFC_I=_get_DIFC_I(mol_h2o, atm_id = i, cartesian=True)
    titulos1_row1=['FC01xx'+"{:02.0f}".format(i+1), 'FC01xy'+"{:02.0f}".format(i+1), 'FC01xz'+"{:02.0f}".format(i+1)]
    titulos1_row2=['FC01yx'+"{:02.0f}".format(i+1), 'FC01yy'+"{:02.0f}".format(i+1), 'FC01yz'+"{:02.0f}".format(i+1)]
    titulos1_row3=['FC01zx'+"{:02.0f}".format(i+1), 'FC01zy'+"{:02.0f}".format(i+1), 'FC01zz'+"{:02.0f}".format(i+1)]

    DIFC_II=_get_DIFC_II(mol_h2o, atm_id = i, cartesian=True)
    titulos2_row1=['FC02xx'+"{:02.0f}".format(i+1), 'FC02xy'+"{:02.0f}".format(i+1), 'FC02xz'+"{:02.0f}".format(i+1)]
    titulos2_row2=['FC02yx'+"{:02.0f}".format(i+1), 'FC02yy'+"{:02.0f}".format(i+1), 'FC02yz'+"{:02.0f}".format(i+1)]
    titulos2_row3=['FC02zx'+"{:02.0f}".format(i+1), 'FC02zy'+"{:02.0f}".format(i+1), 'FC02zz'+"{:02.0f}".format(i+1)]


    with open("shi.bin", "ab") as file:
        file.write(titulos1_row1[0].encode('ascii'))
        numpy.array(DIFC_I[0,0,:,:], dtype=numpy.float64).transpose(1,0).tofile(file)
        file.write(titulos1_row1[1].encode('ascii'))
        numpy.array(DIFC_I[0,1,:,:], dtype=numpy.float64).transpose(1,0).tofile(file)
        file.write(titulos1_row1[2].encode('ascii'))
        numpy.array(DIFC_I[0,2,:,:], dtype=numpy.float64).transpose(1,0).tofile(file)
        file.write(titulos1_row2[0].encode('ascii'))
        numpy.array(DIFC_I[1,0,:,:], dtype=numpy.float64).transpose(1,0).tofile(file)
        file.write(titulos1_row2[1].encode('ascii'))
        numpy.array(DIFC_I[1,1,:,:], dtype=numpy.float64).transpose(1,0).tofile(file)
        file.write(titulos1_row2[2].encode('ascii'))
        numpy.array(DIFC_I[1,2,:,:], dtype=numpy.float64).transpose(1,0).tofile(file)
        file.write(titulos1_row3[0].encode('ascii'))
        numpy.array(DIFC_I[2,0,:,:], dtype=numpy.float64).transpose(1,0).tofile(file)
        file.write(titulos1_row3[1].encode('ascii'))
        numpy.array(DIFC_I[2,1,:,:], dtype=numpy.float64).transpose(1,0).tofile(file)
        file.write(titulos1_row3[2].encode('ascii'))
        numpy.array(DIFC_I[2,2,:,:], dtype=numpy.float64).transpose(1,0).tofile(file)

        file.write(titulos2_row1[0].encode('ascii'))
        numpy.array(DIFC_II[0,0,:,:], dtype=numpy.float64).transpose(1,0).tofile(file)
        file.write(titulos2_row1[1].encode('ascii'))
        numpy.array(DIFC_II[0,1,:,:], dtype=numpy.float64).transpose(1,0).tofile(file)
        file.write(titulos2_row1[2].encode('ascii'))
        numpy.array(DIFC_II[0,2,:,:], dtype=numpy.float64).transpose(1,0).tofile(file)
        file.write(titulos2_row2[0].encode('ascii'))
        numpy.array(DIFC_II[1,0,:,:], dtype=numpy.float64).transpose(1,0).tofile(file)
        file.write(titulos2_row2[1].encode('ascii'))
        numpy.array(DIFC_II[1,1,:,:], dtype=numpy.float64).transpose(1,0).tofile(file)
        file.write(titulos2_row2[2].encode('ascii'))
        numpy.array(DIFC_II[1,2,:,:], dtype=numpy.float64).transpose(1,0).tofile(file)
        file.write(titulos2_row3[0].encode('ascii'))
        numpy.array(DIFC_II[2,0,:,:], dtype=numpy.float64).transpose(1,0).tofile(file)
        file.write(titulos2_row3[1].encode('ascii'))
        numpy.array(DIFC_II[2,1,:,:], dtype=numpy.float64).transpose(1,0).tofile(file)
        file.write(titulos2_row3[2].encode('ascii'))
        numpy.array(DIFC_II[2,2,:,:], dtype=numpy.float64).transpose(1,0).tofile(file)
