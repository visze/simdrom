package de.charite.compbio.simdrom.cli.exception;

/**
 * Exception if an interval (chr1:1212-12121) is formatted wrong.
 * 
 * @author Max Schubach <max.schubach@charite.de>
 *
 */
public class WrongIntervalFormatException extends Exception {

	/**
	 * serial key for serialization
	 */
	private static final long serialVersionUID = 5514299661602182053L;

	/**
	 * Constructor that adds the misspelled string of the interval to the
	 * message.
	 * 
	 * @param misspelled
	 *            String of interval
	 */
	public WrongIntervalFormatException(String input) {
		super(buildMessage(input));

	}

	private static String buildMessage(String input) {
		StringBuilder message = new StringBuilder("Could not parse interval ");
		message.append(input);
		message.append("! Format should be like chr2:121312-121312.");
		return message.toString();
	}

}
