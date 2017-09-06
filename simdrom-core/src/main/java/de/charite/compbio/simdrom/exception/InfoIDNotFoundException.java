/**
 * 
 */
package de.charite.compbio.simdrom.exception;

/**
 * @author max
 *
 */
public class InfoIDNotFoundException extends SIMdromException {

	/**
	 * 
	 */
	private static final long serialVersionUID = -1171419629742433654L;

	/**
	 * 
	 */
	public InfoIDNotFoundException(String id) {
		super("Cannot find INFO ID '" + id
				+ "' in the Head of the VCF file! Every ID have to be present in the VCF header.");
	}

}
