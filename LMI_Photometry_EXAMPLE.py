import LMI_Photometry as LMI



LMI.Data_Reduction(directory='C:\\Users\\caleb\\Documents\\PYTHON\\U_Chicago_Research\\Data\\LMI_2015Jun02\\LMI.20150602',
                   filters={'R':0.1, 'I':0.2},
                   targets={'1RXS_J091744.5+461229' : '2MASS_J09174473+4612246',
                            'PG_1633+099' : 'PG_1633+099'}
                   )

LMI.Aperture_Photometry(directory='C:\\Users\\caleb\\Documents\\PYTHON\\U_Chicago_Research\\Data\\LMI_2015Jun02\\LMI.20150602\\pipeline_out',
                        ap_radius=30,
                        standards={'PG_1633+099' : ['PG_1633+099', 'GSC_00964-01181', 'GSC_00964-01614']}
                        )

LMI.Convert_Magnitudes(directory='C:\\Users\\caleb\\Documents\\PYTHON\\U_Chicago_Research\\Data\\LMI_2015Jun02\\LMI.20150602\\pipeline_out',
                       filters=['I', 'R']
                       )