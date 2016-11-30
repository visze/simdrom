package de.charite.compbio.simdrom.filter;

import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.List;
import java.util.Optional;
import java.util.regex.Pattern;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.CommonInfo;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

/**
 * An implementation to use the VCF-Info column with a key {@link #info} and a value {@link #type} to filter out
 * variants that do not match to key=value.
 * 
 * @author <a href="mailto:max.schubach@charite.de">Max Schubach</a>
 *
 */
public abstract class AInfoFieldFilter extends AFilter {

	/**
	 * key in INFO column
	 */
	private String info;
	/**
	 * Value of the {@link #info} in INFO column
	 */
	private Object type;

	/**
	 * Default constructor
	 * 
	 * @param info
	 *            The key in INFO column
	 * @param type
	 *            Value of the info input in the INFO column
	 */
	public AInfoFieldFilter(String info, Object type) {
		super(FilterType.INFO_FIELD_FILTER);
		this.info = info;
		this.type = type;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.charite.compbio.simdrom.filter.IFilter#filter(htsjdk.variant. variantcontext.VariantContext)
	 */
	@Override
	public Optional<VariantContext> filter(Optional<VariantContext> optional_vc) {
		if (optional_vc.isPresent()) {
			VariantContext vc = optional_vc.get();
			CommonInfo infoField = vc.getCommonInfo();
			if (infoField.hasAttribute(info)) {
				Object val = infoField.getAttribute(info);
				if (val.getClass().isArray())
					return getVariantContextFromArray(vc, val);
				else if (val instanceof List)
					return getVariantContextFromArray(vc, ((List<?>) val).toArray());
				else if (compareInfoType(type, val))
					return optional_vc;
				else {
					return Optional.empty();
				}
			}
		}
		return optional_vc;
	}

	private Optional<VariantContext> getVariantContextFromArray(VariantContext vc, Object infoField) {

		final int length = Array.getLength(infoField);

		// check if we have the same number of attributes than the same number
		// of ALT alleles.
		if (length == vc.getAlternateAlleles().size()) {
			List<Allele> alleles = new ArrayList<Allele>();
			alleles.add(vc.getReference());
			for (int i = 0; i < length; i++) {
				if (compareInfoType(type, Array.get(infoField, i)))
					alleles.add(vc.getAlternateAllele(i));
			}
			// no allele matches, return null
			if (alleles.size() <= 1)
				return Optional.empty();
			else {
				List<Genotype> genotypes = new ArrayList<>();
				for (String sample : vc.getSampleNames()) {
					Genotype g = vc.getGenotype(sample);
					
					List<Allele> g_alleles = new ArrayList<>();
					for (Allele allele : g.getAlleles()) {
						if (alleles.contains(allele))
							g_alleles.add(allele);
						else {
							g_alleles.add(Allele.NO_CALL);
						}
					}
					GenotypeBuilder gb = new GenotypeBuilder().phased(g.isPhased()).alleles(g_alleles).name(sample);
					
					genotypes.add(gb.make());
				}
				
				return Optional.of(new VariantContextBuilder(vc).alleles(alleles).genotypes(genotypes).make());
			}
		} else { // hack, now we add if one allele matches it.
			for (int i = 0; i < length; i++) {
				if (compareInfoType(type, Array.get(infoField, i)))
					return Optional.of(vc);
			}
		}
		return Optional.empty();

	}

	protected abstract boolean compareInfoType(Object should, Object is);

	protected int compare(Object should, Object is) {
		if (should instanceof String) {
			return ((String) should).compareTo((String) is);
		} else if (should instanceof Integer)
			return ((Integer) should).compareTo(Integer.valueOf((String) is));
		else if (should instanceof Double)
			return ((Double) should).compareTo(Double.valueOf((String) is));

		throw new RuntimeException("Cannot compare object " + should + " with object " + is);
	}

	private static boolean isDouble(String string) {
		return Pattern.matches("^\\d+(.\\d+)?$", string);
	}

	private static boolean isInt(String string) {
		return Pattern.matches("^\\d+$", string);
	}

}
