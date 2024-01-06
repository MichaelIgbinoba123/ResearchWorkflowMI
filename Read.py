import netCDF4 as nc

# Specify the path to the netCDF file
file_path = r'D:\Noru2022.nc'

# Open the netCDF file
with nc.Dataset(file_path, 'r') as dataset:
    # Display global attributes
    print("Global attributes:")
    for attr_name in dataset.ncattrs():
        print(f"{attr_name}: {getattr(dataset, attr_name)}")

    # Display variables and their attributes
    print("\nVariables and their attributes:")
    for var_name, variable in dataset.variables.items():
        print(f"\nVariable: {var_name}")
        print(f"Shape: {variable.shape}")
        print("Attributes:")
        for attr_name in variable.ncattrs():
            print(f"  {attr_name}: {getattr(variable, attr_name)}")