package de.charite.compbio.simdrom.filter;

/**
 * Enum of different Filters ({@link IFilter}) used on a VCF-File
 * 
 * @author <a href="mailto:max.schubach@charite.de">Max Schubach</a>
 *
 */
public enum FilterType {
	/**
	 * Filter that uses the VCF-INFO field with a specific key=value tag to remove variants without this key=value.
	 */
	INFO_FIELD_FILTER,
	/**
	 * Filter the clinVar VCF file looking for pathogenic variants. It takes the Clinical allele , the clinical
	 * significance, the clincial origin and the clinical db/in into account.
	 */
	CLINVAR_FILTER,
	/**
	 * Removes all gonosomes. Sometimes Gonosomes might be not needed for some applications.   
	 */
	GONOSOME_FILTER;
}
