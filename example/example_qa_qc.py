from GEOAkaze import qa_geoakaze

geoakaze_nc_fld = "/n/holylfs04/LABS/wofsy_lab/Lab/MethaneAIR/level1/RF06_V3/GEOAkazeNC_O2"
geoakaze_kmz_fld = "/n/holylfs04/LABS/wofsy_lab/Lab/MethaneAIR/level1/RF06_V3/GEOAkazeKMZ_O2"
L1b_folder = "/n/holylfs04/LABS/wofsy_lab/Lab/MethaneAIR/level1/RF06_V3/O2_NATIVE"
output_pdf = "./report_ortho_RF08.pdf"
temp_fld = "./temp_qa_qc/"

qa_obj = qa_geoakaze(geoakaze_nc_fld,geoakaze_kmz_fld,L1b_folder,output_pdf,temp_fld)

#qa_obj.gray_scale(random_selection_n=10)
#qa_obj.histogram()
#qa_obj.trajectory()
#qa_obj.counter()

qa_obj.topdf()
