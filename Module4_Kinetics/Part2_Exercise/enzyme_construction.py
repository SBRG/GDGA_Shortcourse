"""Contains functions for facilitating the E.coli glycolytic case study."""

import os
import re
from operator import attrgetter, iconcat

from six import iteritems

import sympy as sym

from mass import MassMetabolite
from mass.enzyme_modules import (
    EnzymeModule, EnzymeModuleForm, EnzymeModuleReaction)
from mass.util.expressions import strip_time

# Pre-compiled regex
compartment_finder_re = re.compile(r"\s*\[([A-Za-z])\]")
rateconst_re = re.compile(r"k_(\S*)_(fwd|rev)\Z")
bound_metabolites_re = re.compile(r"&|@|#")


def format_percent_str(percent):
    return str(int(round(percent * 100, 0))).replace(".", "")


def prefix_number_id(id_str):
    """Prefix identifiers that start with numbers."""
    if re.match(r"^\d", id_str):
        id_str = "_" + id_str
    return id_str


def make_path(*args):
    """Combine path arguments and return absolute filepath."""
    return os.path.abspath(os.path.join(*args))

def fix_equation_variables(enzyme_module, equation_str):
    """Make corrections for equation strings in order to sympify."""
    # Correct compartments paranthesis to brackets
    for cid in enzyme_module.compartments:
        c_str = r"\({0}\)".format(cid)
        equation_str = re.sub(
            c_str, lambda match: "[{0}]".format(cid), equation_str)

        vol_str = "_".join(("param", "Volume", cid))
        if vol_str in equation_str:
            equation_str = re.sub(vol_str, "1", equation_str)

        e_tot_str = "_".join(("param", enzyme_module.id, "total"))
        if e_tot_str in equation_str:
            equation_str = re.sub(
                e_tot_str, enzyme_module.enzyme_total_symbol_str, equation_str)

    for dtype in ["species", "kf", "kr"]:
        for old_id, new_id in iteritems(enzyme_module.id_map[dtype]):
            if dtype == "species":
                old_id = old_id.replace("[", r"\[").replace("]", r"\]")
            equation_str = re.sub(old_id, lambda match: new_id, equation_str)
    return equation_str


def calculate_enzyme_total(enzyme_module, concentrations, steady_state_flux):
    """Calculate the total amount of enzyme necessary to sustain given flux."""
    rate_eq = sym.Eq(
        enzyme_module.enzyme_rate_equation.subs(concentrations),
        steady_state_flux)
    e_total = list(sym.solveset(
        rate_eq.subs(concentrations), enzyme_module.enzyme_total_symbol_str))
    e_total = float(e_total.pop())
    return e_total


def make_enzyme_module_from_dir(enzyme_id, steady_state_flux=None,
                                metabolite_concentrations=None,
                                path_to_dir=None, **kwargs):
    """Create an enzyme module from a directory of text files."""
    enzyme_module = EnzymeModule(enzyme_id)
    # Add temporary attribute to track ID mapping changes
    enzyme_module.id_map = {}
    # Set flux through enzyme
    enzyme_module.enzyme_rate = steady_state_flux

    def path_to_enzyme_file(type_str):
        """Return filepath to an enzyme data file."""
        enzyme_dir = os.path.join(path_to_dir, enzyme_module.id)
        filename = "_".join((type_str, enzyme_module.id.lower()))
        if filename + ".txt" in os.listdir(enzyme_dir):
            filename += ".txt"
        return os.path.join(enzyme_dir, filename)

    # Make species and add to model
    metabolites_from_txt(enzyme_module, path_to_enzyme_file("species"))

    # Make reactions and add to model
    reactions_from_txt(enzyme_module, path_to_enzyme_file("reactions"))

    # Set GPR if given.
    if kwargs.get("enzyme_gpr", None) is not None:
        for reaction in enzyme_module.enzyme_module_reactions:
            reaction.gene_reaction_rule = kwargs.get("enzyme_gpr", None)
    # Determine rate constant values and add to model
    rateconsts_from_txt(enzyme_module, kwargs.get("kcluster", 1),
                        path_to_enzyme_file("rateconst_labels"),
                        path_to_enzyme_file("rateconst_clusters"))
    # Set rate law equation for flux through enzyme.
    kinetic_rate_law_from_txt(enzyme_module, path_to_enzyme_file("rateLaw"))
    # Calculate the expected total amount of enzyme based on the enzyme flux
    metabolite_concentrations = {
        m.id: ic for m, ic in iteritems(metabolite_concentrations)
        if m.id in enzyme_module.metabolites}
    enzyme_module.enzyme_concentration_total = calculate_enzyme_total(
        enzyme_module, metabolite_concentrations, steady_state_flux)

    # Get steady state concentrations for enzyme module forms
    for old_specie, specie in iteritems(enzyme_module.id_map["species"]):
        try:
            specie = enzyme_module.enzyme_module_forms.get_by_id(specie)
        except KeyError:
            specie = enzyme_module.metabolites.get_by_id(specie)
            specie.ic = metabolite_concentrations[specie.id]
        else:
            steady_state_concentrations_from_txt(
                enzyme_module, path_to_enzyme_file("equation_" + old_specie),
                specie, metabolite_concentrations,
                zero_tol=kwargs.get("zero_tol", 1e-15))
            specie._repair_bound_obj_pointers()
            specie.generate_form_formula(update_enzyme=True)
            specie.generate_form_charge(update_enzyme=True)

    # Remove ID map before returning module
    del enzyme_module.id_map

    return enzyme_module


def metabolites_from_txt(enzyme_module, species_filepath):
    """Create metabolites from a text file."""
    # Create ID map dict for species
    enzyme_module.id_map["species"] = {}

    # Read metabolites from files
    with open(species_filepath, "r") as f:
        original_ids = [l.strip() for l in f.readlines()]

    species = []
    for orig_id in original_ids:
        # Initialize dict for object initialization kwargs and variable for ID
        new_id = orig_id
        # Find compartment, if any.
        compartment = compartment_finder_re.findall(orig_id)
        if compartment is not None:
            # Correct isomer in ID
            if "D" in compartment or "L" in compartment:
                isomer = compartment.pop(
                    compartment.index("D" if "D" in compartment else "L"))
                new_id = re.sub(r"\[" + isomer + r"\]", "_" + isomer, new_id)

            # One compartment remaining indicates compartment for specie.
            if len(compartment) == 1:
                compartment = compartment.pop()
                # Correct compartment in ID
                new_id = compartment_finder_re.sub("_" + compartment, new_id)

        # Make ID corrections
        new_id = prefix_number_id(new_id)
        # Create specie
        if new_id.startswith("E_" + enzyme_module.id):
            bound_metabolites = {}
            if bound_metabolites_re.search(new_id):
                # Correct IDs on bound species and store number bound
                for bound_specie_id in bound_metabolites_re.split(new_id)[1:]:
                    if compartment is not None:
                        bound_specie_id += "_" + compartment
                    bound_specie_id = prefix_number_id(bound_specie_id)
                    if bound_specie_id not in bound_metabolites:
                        bound_metabolites[bound_specie_id] = 0
                    bound_metabolites[bound_specie_id] += 1
                # Correct ID for enzyme form.
                new_id = bound_metabolites_re.sub("_", new_id)

            specie = EnzymeModuleForm(
                id_or_specie=new_id,
                enzyme_module_id=enzyme_module.id,
                compartment=compartment)
            specie._bound_metabolites = bound_metabolites
        else:
            specie = MassMetabolite(
                id_or_specie=new_id,
                compartment=compartment)
        species.append(specie)
        # Store ID mapping
        enzyme_module.id_map["species"][orig_id] = new_id
    # Add MassMetabolites and EnzymeModuleForms to enzyme module
    enzyme_module.add_metabolites(species)

    return None


def reactions_from_txt(enzyme_module, reactions_filepath):
    """Create reactions from a text file."""
    # Create ID map dict for reactions
    enzyme_module.id_map["reactions"] = {}

    # Read reactions from files
    with open(reactions_filepath, "r") as f:
        reaction_strings = [l.strip() for l in f.readlines()]

    for reaction_str in reaction_strings:
        # Split reaction ID from reaction formula
        orig_id, reaction_str = reaction_str.split(r": ")
        new_id = prefix_number_id(re.sub(r"\$", "_", orig_id))
        # Make EnzymeModuleReaction object and add to model
        reaction = EnzymeModuleReaction(
            id_or_reaction=new_id,
            enzyme_module_id=enzyme_module.id)
        enzyme_module.add_reactions([reaction])

        # Fix specie IDs in reaction string
        species_list = iconcat(*[
            s.split(r" + ") for s in reaction_str.strip().split(r" <=> ")])
        for specie_id in species_list:
            reaction_str = re.sub(
                specie_id.replace(r"[", r"\[").replace(r"]", r"\]"),
                enzyme_module.id_map["species"][specie_id],
                reaction_str, 1)
        # Build reaction from string.
        reaction.build_reaction_from_string(reaction_str)
        # Store ID mapping
        enzyme_module.id_map["reactions"][orig_id] = new_id

    return None


def rateconsts_from_txt(enzyme_module, kcluster, labels_filepath,
                        clusters_filepath):
    """Get values of forward and reverse rate constants from text files."""
    enzyme_module.id_map["kf"] = {}
    enzyme_module.id_map["kr"] = {}

    # Read rate constant labels from files
    with open(labels_filepath, "r") as f:
        rateconst_labels = [l.strip() for l in f.readlines()]
    # Read rate constant values from files
    with open(clusters_filepath, "r") as f:
        lines = f.readlines()
    try:
        values_txt = lines[kcluster - 1]
    except IndexError:
        raise IndexError(
            "Invalid kcluster. Must be an int between [1, {}], ".format(
                int(len(lines))))
    values_txt = values_txt.strip().lstrip("{").rstrip("}")

    # Assign values to rate constants
    to_unify = {}
    rateconst_values = {}
    for rateconst, value in zip(rateconst_labels, values_txt.split(", ")):
        match = rateconst_re.search(rateconst)
        rid, direction = match.groups()
        key = {"fwd": "kf", "rev": "kr"}[direction]
        try:
            rxn = enzyme_module.reactions.get_by_id(
                enzyme_module.id_map["reactions"][rid])
        except KeyError:
            new_rateconst = "_".join((key, rid))
            to_unify[rid] = [r for r in enzyme_module.reactions if rid in r.id]
        else:
            new_rateconst = getattr(rxn, "_".join((key, "str")))

        try:
            value = float(value)
        except ValueError:
            value = float(value.replace("*^", "e"))
        enzyme_module.id_map[key][rateconst] = new_rateconst
        rateconst_values[new_rateconst] = value * 3600  # s -> hrs

    # Handle any symmetry model rate laws
    for rid, rxn_list in iteritems(to_unify):
        # Unity rate parameters
        enzyme_module.unify_rate_parameters(rxn_list, rid, rate_type=2)
        # Set symmetry coefficients for rate constants in rates
        kf_str, kr_str = ["_".join((k, rid)) for k in ["kf", "kr"]]
        for i, rxn in enumerate(sorted(rxn_list, key=attrgetter("id"))):
            rate_equation = str(strip_time(rxn.rate))
            rate_equation = rate_equation.replace(
                kr_str, "{0}*{1}".format(str(i + 1), kr_str))
            rate_equation = rate_equation.replace(
                kf_str, "{0}*{1}".format(str(len(rxn_list) - i), kf_str))
            # Add as a custom rate
            enzyme_module.add_custom_rate(
                rxn, rate_equation, custom_parameters={
                    kf_str: None, kr_str: None})

    # Update enzyme module with values
    enzyme_module.update_parameters(rateconst_values)

    return None


def kinetic_rate_law_from_txt(enzyme_module, rate_law_filepath):
    """Get the symbolic rate law equation from a text file."""
    # Read rate constant labels from files
    with open(rate_law_filepath, "r") as f:
        rate_law = f.read()
    enzyme_module.get_rate_expressions(rate_type=2, update_reactions=True)
    rate_law = fix_equation_variables(enzyme_module, rate_law)
    enzyme_module.enzyme_rate_equation = sym.sympify(
        rate_law, locals=enzyme_module._get_all_parameters())

    return None


def steady_state_concentrations_from_txt(enzyme_module, eq_filepath, specie,
                                         concentrations, zero_tol):
    """Get an equation for enzyme form concentration from a text file."""
    # Read rate constant labels from files
    with open(eq_filepath, "r") as f:
        ss_eq = f.read()
    ss_eq = fix_equation_variables(enzyme_module, ss_eq)
    ss_eq = sym.sympify(ss_eq, locals=enzyme_module._get_all_parameters())
    e_tot = {
        enzyme_module.id + "_Total": enzyme_module.enzyme_concentration_total}
    ic = float(ss_eq.subs(concentrations).subs(e_tot))
    if abs(ic) <= zero_tol:
        ic = 0
    specie.ic = ic

    return None