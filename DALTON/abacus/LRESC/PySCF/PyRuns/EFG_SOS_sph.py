#import os 
#os.system("export PYTHONPATH=/home/juanjoaucar/Archivos_Juan/pyscf_new:$PYTHONPATH")
#import pyscf
from pyscf import gto,scf
import numpy

def _get_new_integrals(mol, atm_id, zeta = 9999999999, cartesian = False):
#    nao = mol.nao
    nao_cart = mol.intor('int1e_ovlp_cart').shape[0]
    nao_sph = mol.intor('int1e_ovlp_sph').shape[0]
    with mol.with_rinv_zeta(zeta), mol.with_rinv_origin((mol.atom_coord(atm_id))):
#        c = mol.cart2sph_coeff()
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

        #Chequeado para el elemento zzzz (2,2,2,2,:,:) con funciones s y con funciones p:
        core = -ipipiprinvip.transpose(3,1,2,0,5,4)-ipiprinvipip.transpose(0,3,1,2,4,5)-ipipiprinvip.transpose(1,3,2,0,5,4)-ipiprinvipip.transpose(3,0,1,2,4,5)-ipipiprinvip.transpose(1,2,3,0,5,4)-ipipipiprinv.transpose(1,2,3,0,5,4)-ipipiprinvip.transpose(0,1,2,3,4,5)-ipiprinvipip.transpose(0,1,3,2,4,5)

        #Transpose (gives same result)
#        core = -ipipiprinvip.transpose(3,1,0,2,5,4)-ipiprinvipip.transpose(1,2,0,3,4,5)-ipipiprinvip.transpose(1,3,0,2,5,4)-ipiprinvipip.transpose(2,1,0,3,4,5)-ipipiprinvip.transpose(1,0,3,2,5,4)-ipipipiprinv.transpose(2,1,0,3,5,4)-ipipiprinvip.transpose(2,1,0,3,4,5)-ipiprinvipip.transpose(1,0,2,3,4,5)


    return core

def _get_GiZZj(mol,atm_id,i,j, zeta = 9999999999, cartesian = False):
    #Gets the transformed matrix
    matriz_xx=_get_new_integrals(mol,atm_id,zeta,cartesian)[i,0,0,j,:,:]
    matriz_yy=_get_new_integrals(mol,atm_id,zeta,cartesian)[i,1,1,j,:,:]
    matriz_zz=_get_new_integrals(mol,atm_id,zeta,cartesian)[i,2,2,j,:,:]
    return 1/3*(2*matriz_zz-matriz_xx-matriz_yy)



#base = "dyall_cv3z"
base = "dyall_cv2z"

mol_h2o = gto.M(atom="molecula.xyz",basis = base)
mol_h2o, ctr_coeff = mol_h2o.to_uncontracted_cartesian_basis()

#orbitales = numpy.array([mol_h2o.nao])
naos_sph = mol_h2o.intor('int1e_ovlp_sph').shape[0]
nao_cart = mol_h2o.intor('int1e_ovlp_cart').shape[0]
nro_operadores=numpy.array([3*mol_h2o.natm])
#print("naos_cart:",orbitales)
print("naos_cart:",nao_cart)
print("naos_sph:",naos_sph)
print("nro de operadores a grabar: ",nro_operadores)
cartesian = False
if (cartesian):
   print("Cartesian basis set used at PySCF")
else:
   print("Spherical basis set used at PySCF")

with open("efg.bin", "wb") as file:
    numpy.array(nro_operadores[0], dtype=numpy.int32).tofile(file)
    if (cartesian):
       numpy.array(nao_cart, dtype=numpy.int32).tofile(file)
    else:
       numpy.array(naos_sph, dtype=numpy.int32).tofile(file)


for i in range(mol_h2o.natm):
    titulos=['newopx'+"{:02.0f}".format(i+1), 'newopy'+"{:02.0f}".format(i+1), 'newopz'+"{:02.0f}".format(i+1)]
    print(titulos)
    valorzeta=EXPONENTEMANUAL #Para variar en script
#    valorzeta=10E+11 #Point Model
#    valorzeta=1E+05


#   Minus sign added due to p
    nuevo_x=-(_get_GiZZj(mol_h2o,i,1,2, zeta=valorzeta, cartesian = cartesian)-_get_GiZZj(mol_h2o,i,2,1, zeta=valorzeta, cartesian = cartesian))
    nuevo_y=-(_get_GiZZj(mol_h2o,i,2,0, zeta=valorzeta, cartesian = cartesian)-_get_GiZZj(mol_h2o,i,0,2, zeta=valorzeta, cartesian = cartesian))
    nuevo_z=-(_get_GiZZj(mol_h2o,i,0,1, zeta=valorzeta, cartesian = cartesian)-_get_GiZZj(mol_h2o,i,1,0, zeta=valorzeta, cartesian = cartesian))

    with open("efg.bin", "ab") as file:
        file.write(titulos[0].encode('ascii'))
        numpy.array(nuevo_x, dtype=numpy.float64).transpose(1,0).tofile(file)
        print("grabando en efg.bin el operador ", titulos[0])
        file.write(titulos[1].encode('ascii'))
        numpy.array(nuevo_y, dtype=numpy.float64).transpose(1,0).tofile(file)
        file.write(titulos[2].encode('ascii'))
        numpy.array(nuevo_z, dtype=numpy.float64).transpose(1,0).tofile(file)
