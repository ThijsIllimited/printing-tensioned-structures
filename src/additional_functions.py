def replace_brackets(line, temperature_settings):
    """"Replace the brackets and the brackets contents with the corresponding variable from the dictionary"""
    for key, value in temperature_settings.items():
        key_b = '[' + key + ']'
        if key_b in line:
            line = line.replace(key_b, str(value))
    return line