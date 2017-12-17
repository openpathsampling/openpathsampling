
# TODO: maybe move the class initialization and the variable initialization
# into different levels to allow only class initialization to be logged?
def initialization_logging(logger, obj, entries):
    # this works with either a list of attributes of a dictionary of
    # variable name : attribute value
    try:
        entries.keys()
        working_dict = entries
    except AttributeError:
        working_dict = {}
        for entry in entries:
            working_dict[entry] = obj.__dict__[entry]

    logger.info("Initializing <%s> (%s)",
                str(hex(id(obj))), str(obj.__class__.__name__))
    for entry in working_dict.keys():
        logger.info("Parameter: %s : %s", str(entry), str(working_dict[entry]))
